#ifndef INC_LRENDER_BINBVH_H__
#define INC_LRENDER_BINBVH_H__
/**
@file BinBVH.h
@author t-sakai
@date 2015/09/09 create
*/
#include "lrender.h"
#include <fstream>

namespace lrender
{
    template<class PrimitiveType, class PrimitivePolicy = lrender::PrimitivePolicy<PrimitiveType> >
    class BinBVH
    {
    public:
        static const f32 Epsilon;
        static const s32 MinLeafPrimitives = 15;
        static const s32 NumBins = 32;
        static const s32 MaxBinningDepth = 16;

        struct Node
        {
            static const u32 LeafMask = ~(((u32)-1)>>1);
            static const u32 EmptyMask = LeafMask;

            static const s32 MaxNumLeafPrimitiveShift = 4;
            static const s32 LeafPrimitiveMask = (0x01U<<MaxNumLeafPrimitiveShift)-1;
            static const u32 MaxNumPrimitives = (~(LeafMask|EmptyMask))>>MaxNumLeafPrimitiveShift;

            bool isLeaf() const
            {
                return (LeafMask&child_) != 0;
            }

            bool isEmpty() const
            {
                return (EmptyMask == child_);
            }

            void setLeaf(u32 primitiveIndex, u32 numPrimitives)
            {
                LASSERT(0<numPrimitives);
                child_ =  LeafMask | ((primitiveIndex&MaxNumPrimitives) << MaxNumLeafPrimitiveShift) | (numPrimitives&LeafPrimitiveMask);
            }

            u32 getPrimitiveIndex() const
            {
                return (child_ >> MaxNumLeafPrimitiveShift) & MaxNumPrimitives;
            }

            u32 getNumPrimitives() const
            {
                return (child_&0x0FU);
            }

            BBox bbox_;
            u32 child_;
        };

        BinBVH();
        ~BinBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);

        s32 getDepth() const{ return depth_;}
        void print(const char* filename);
    private:
        BinBVH(const BinBVH&);
        BinBVH& operator=(const BinBVH&);

        inline void getBBox(BBox& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 start, s32 numPrimitives, s32 nodeIndex, s32 depth);
        void split(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 invArea, s32 start, s32 numPrimitives, s32 nodeIndex);
        void splitMid(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, s32 start, s32 numPrimitives, s32 nodeIndex);
        void splitBinned(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 area, s32 start, s32 numPrimitives, s32 nodeIndex);

        f32 SAH_KI_;
        f32 SAH_KT_;
        const PrimitiveType* primitives_;

        s32 depth_;
        vector_arena<Node> nodes_;
        vector_arena<s32> primitiveIndices_;
        vector_arena<f32> primitiveCentroids_;
        vector_arena<BBox> primitiveBBoxes_;
        s32* stack_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    const f32 BinBVH<PrimitiveType, PrimitivePolicy>::Epsilon = 1.0e-6f;

    template<class PrimitiveType, class PrimitivePolicy>
    BinBVH<PrimitiveType, PrimitivePolicy>::BinBVH()
        :SAH_KI_(1.5f)
        ,SAH_KT_(1.0f)
        ,primitives_(NULL)
        ,depth_(0)
        ,stack_(NULL)
    {
        nodes_.resize(1);
        nodes_[0].bbox_.zero();
        nodes_[0].child_ = Node::EmptyMask;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    BinBVH<PrimitiveType, PrimitivePolicy>::~BinBVH()
    {
        LIME_FREE(stack_);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
    {
        f32 depth = logf(static_cast<f32>(numPrimitives) / MinLeafPrimitives) / logf(2.0f);
        s32 numNodes = static_cast<s32>(powf(2.0f, depth) + 0.5f);
        nodes_.reserve(numNodes);
        nodes_.resize(1);
        primitiveIndices_.resize(numPrimitives);
        primitiveCentroids_.resize(numPrimitives*3);
        primitiveBBoxes_.resize(numPrimitives);

        primitives_ = primitives;
        nodes_[0].bbox_.setInvalid();

        //各primitiveのcentroid, bboxを事前計算
        f32* centroidX = &primitiveCentroids_[0];
        f32* centroidY = centroidX + numPrimitives;
        f32* centroidZ = centroidY + numPrimitives;
        for(s32 i=0; i<numPrimitives; ++i){
            primitiveIndices_[i] = i;

            Vector3 centroid = PrimitivePolicy::getCentroid(primitives_[i]);
            centroidX[i] = centroid.x_;
            centroidY[i] = centroid.y_;
            centroidZ[i] = centroid.z_;
            primitiveBBoxes_[i] = PrimitivePolicy::getBBox(primitives_[i]);
            nodes_[0].bbox_.extend( primitiveBBoxes_[i] );
        }

        s32 prevDepth = depth_;
        depth_ = 0;
        recursiveConstruct(0, numPrimitives, 0, 1);

        primitiveCentroids_.clear();
        primitiveBBoxes_.clear();

        //intersectのためにスタック準備
        if(prevDepth<depth_){
            LIME_FREE(stack_);
            stack_ = (s32*)LIME_MALLOC(sizeof(s32)*depth_);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    inline void BinBVH<PrimitiveType, PrimitivePolicy>::getBBox(BBox& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[primitiveIndices_[i]]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 start, s32 numPrimitives, s32 nodeIndex, s32 depth)
    {
        if (0 == numPrimitives){
            nodes_[nodeIndex].child_ = Node::EmptyMask;
            depth_ = maximum(depth, depth_);
            return;
        }
        if (numPrimitives <= MinLeafPrimitives){
            nodes_[nodeIndex].setLeaf(start, numPrimitives);
            depth_ = maximum(depth, depth_);
            return;
        }

        s32 num_l, num_r;
        BBox bbox_l, bbox_r;
        f32 area = nodes_[nodeIndex].bbox_.halfArea();
        s32 axis = 0;

        if(MaxBinningDepth<depth || area<=Epsilon){
            splitMid(axis, num_l, num_r, bbox_l, bbox_r, start, numPrimitives, nodeIndex);

        }else if(numPrimitives<NumBins){
            split(axis, num_l, num_r, bbox_l, bbox_r, 1.0f/area, start, numPrimitives, nodeIndex);

        } else{
            splitBinned(axis, num_l, num_r, bbox_l, bbox_r, area, start, numPrimitives, nodeIndex);
        }

        s32 childIndex = nodes_.size();
        nodes_.push_back(Node());
        nodes_.push_back(Node());

        nodes_[nodeIndex].child_ = childIndex;

        nodes_[childIndex].bbox_ = bbox_l;
        nodes_[childIndex+1].bbox_ = bbox_r;

        ++depth;
        recursiveConstruct(start, num_l, nodes_[nodeIndex].child_, depth);
        recursiveConstruct(start+num_l, num_r, nodes_[nodeIndex].child_+1, depth);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HitRecord BinBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
    {
        HitRecord hitRecord;
        hitRecord.t_ = ray.t_;
        hitRecord.primitive_ = NULL;
        s32 top = 0;

        s32 currentNode = 0;
        for(;;){
            Node& node = nodes_[currentNode];

            //f32 tmin, tmax;
            //bool hit = node.bbox_.testRay(tmin, tmax, ray);
            //if(hit){
            //    printf("%f, %f, %f\n", ray.direction_.x_, ray.direction_.y_, ray.direction_.z_);
            //}
            if(node.isLeaf()){
                s32 index = node.getPrimitiveIndex();
                s32 num = node.getNumPrimitives();

                for(s32 i=0; i<num; ++i){
                    f32 t;
                    s32 idx = primitiveIndices_[index+i];
                    if(!primitives_[idx].testRay(t, ray)){
                        continue;
                    }
                    if(HitEpsilon < t && t < hitRecord.t_){
                        ray.t_ = t;
                        hitRecord.t_ = t;
                        hitRecord.primitive_ = &primitives_[idx];
                    }
                }
                if(top <= 0){
                    break;
                }
                currentNode = stack_[--top];

            }else{
                Node& nodeLeft = nodes_[node.child_];
                Node& nodeRight = nodes_[node.child_+1];

                f32 leftMin, leftMax;
                bool hitLeft = nodeLeft.bbox_.testRay(leftMin, leftMax, ray);
                leftMin = maximum(HitEpsilon, leftMin);
                leftMax = minimum(hitRecord.t_, leftMax);

                f32 rightMin, rightMax;
                bool hitRight =  nodeRight.bbox_.testRay(rightMin, rightMax, ray);
                rightMin = maximum(HitEpsilon, rightMin);
                rightMax = minimum(hitRecord.t_, rightMax);

                if(hitLeft && !hitRight){
                    currentNode = node.child_;
                }else if(!hitLeft && hitRight){
                    currentNode = node.child_ + 1;
                }else if(hitLeft && hitRight){
                    currentNode = node.child_;
                    s32 farNode = node.child_+1;
                    if(rightMin<leftMin){
                        swap(currentNode, farNode);
                    }
                    stack_[top++] = farNode;
                }else{
                    if(top<=0){
                        break;
                    }
                    currentNode = stack_[--top];
                }
            }//if(node.isLeaf()){
        }//for(;;){
        return hitRecord;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::split(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 invArea, s32 start, s32 numPrimitives, s32 nodeIndex)
    {
        s32 end = start + numPrimitives;
        s32 mid=start+(numPrimitives >> 1);

        f32 area_l, area_r;
        f32 bestCost = std::numeric_limits<f32>::max();

#if 1
        //SAH, 全ての分割を試す
        axis = nodes_[nodeIndex].bbox_.maxExtentAxis();
        f32* bestCentroids = &primitiveCentroids_[0] + axis*primitiveIndices_.size();
        PrimitivePolicy::insertionsort(numPrimitives, &primitiveIndices_[start], bestCentroids);

        BBox bl, br;
        for(s32 m=start+1; m<end; ++m){
            getBBox(bl, start, m);
            getBBox(br, m, end);

            area_l = bl.halfArea();
            area_r = br.halfArea();
            num_l = m-start;
            num_r = numPrimitives - num_l;

            f32 cost = SAH_KT_ + SAH_KI_*invArea*(area_l*num_l + area_r*num_r);
            if(cost<bestCost){
                mid = m;
                bestCost = cost;
                bbox_l = bl;
                bbox_r = br;
            }
        }

#else
        //SAH, 軸と全ての分割を試す
        axis = 0;
        f32* centroids = &primitiveCentroids_[0];
        f32* bestCentroids = centroids;
        for(s32 curAxis=0; curAxis<3; ++curAxis){
            PrimitivePolicy::insertionsort(numPrimitives, &primitiveIndices_[start], centroids);

            BBox bl, br;
            for(s32 m=start+1; m<end; ++m){
                getBBox(bl, start, m);
                getBBox(br, m, end);

                area_l = bl.halfArea();
                area_r = br.halfArea();
                num_l = m-start;
                num_r = numPrimitives - num_l;

                f32 cost = SAH_KT_ + SAH_KI_*invArea*(area_l*num_l + area_r*num_r);
                if(cost<bestCost){
                    mid = m;
                    bestCost = cost;
                    axis = curAxis;
                    bbox_l = bl;
                    bbox_r = br;
                    bestCentroids = centroids;
                }
            }

            centroids += primitiveIndices_.size();
        }//for(s32 curAxis=0;
        PrimitivePolicy::insertionsort(numPrimitives, &primitiveIndices_[start], bestCentroids);
#endif
        num_l = mid-start;
        num_r = numPrimitives - num_l;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::splitMid(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, s32 start, s32 numPrimitives, s32 nodeIndex)
    {
        //最大の軸を半分に分割
        Node& node = nodes_[nodeIndex];
        axis = node.bbox_.maxExtentAxis();

        s32 end = start + numPrimitives;
        num_l = (numPrimitives >> 1);
        num_r = numPrimitives - num_l;
        s32 mid=start+num_l;

        f32* centroids = &primitiveCentroids_[0] + axis * primitiveIndices_.size();
        PrimitivePolicy::sort(numPrimitives, &primitiveIndices_[start], centroids);

        getBBox(bbox_l, start, mid);
        getBBox(bbox_r, mid, end);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::splitBinned(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 area, s32 start, s32 numPrimitives, s32 nodeIndex)
    {
        LIME_ALIGN16 s32 minBins[NumBins];
        LIME_ALIGN16 s32 maxBins[NumBins];

        __m128 zero = _mm_setzero_ps();

        f32 invArea = 1.0f/area;
        axis = 0;
        s32 end = start + numPrimitives;

        f32* centroids = &primitiveCentroids_[0];

        f32 bestCost = std::numeric_limits<f32>::max();
        s32 midBin = NumBins/2;
        s32 step = static_cast<s32>(::log10f(static_cast<f32>(numPrimitives)));
        Node& node = nodes_[nodeIndex];

        Vector3 extent =node.bbox_.extent();
        Vector3 unit = extent * (1.0f/NumBins);
        for(s32 curAxis=0; curAxis<3; ++curAxis){
            for(s32 i=0; i<NumBins; i+=4){
                _mm_store_ps(reinterpret_cast<f32*>(&minBins[i]), zero);
                _mm_store_ps(reinterpret_cast<f32*>(&maxBins[i]), zero);
            }
            PrimitivePolicy::sort(numPrimitives, &primitiveIndices_[start], centroids);

            f32 invUnit = (absolute(unit[curAxis])<Epsilon)? 0.0f : 1.0f/unit[curAxis];
            f32 bmin = node.bbox_.bmin_[curAxis];

            for(s32 i = start; i < end; i+=step){
                s32 index = primitiveIndices_[i];
                s32 minIndex = minimum(static_cast<s32>(invUnit * (primitiveBBoxes_[index].bmin_[curAxis] - bmin)), NumBins-1);
                s32 maxIndex = minimum(static_cast<s32>(invUnit * (primitiveBBoxes_[index].bmax_[curAxis] - bmin)), NumBins-1);
                ++minBins[minIndex];
                ++maxBins[maxIndex];
            }

            Vector3 e = extent; e[curAxis] = unit[curAxis];
            f32 unitArea = e.halfArea();

            s32 binLeft = 0;
            s32 binRight = NumBins - 1;

            while(minBins[binLeft]<=0){++binLeft;}
            while(maxBins[binRight]<=0){--binRight;}

            s32 n_l = minBins[0];
            s32 n_r = maxBins[0];
            for(s32 i=1; i<=binRight; ++i){
                n_r += maxBins[i];
            }
            for(s32 m=binLeft; m<=binRight; ++m){
                f32 area_l = m*unitArea;
                f32 area_r = (NumBins-m)*unitArea;
                f32 cost = SAH_KT_ + SAH_KI_*invArea*(area_l*n_l + area_r*n_r);
                if(cost<bestCost){
                    midBin = m;
                    bestCost = cost;
                    axis = curAxis;
                }

                n_l += minBins[m];
                n_r -= maxBins[m];
            }

            centroids += primitiveIndices_.size();
        }//for(s32 curAxis=0;

        f32 separate = unit[axis] * (midBin+1) + node.bbox_.bmin_[axis];
        s32 mid = start+(numPrimitives >> 1);
#if 1
        s32 left = start;
        s32 right = end - 1;
        for(;;){
            while(primitiveBBoxes_[ primitiveIndices_[left] ].bmin_[axis] < separate){
                ++left;
            }
            while(separate <= primitiveBBoxes_[ primitiveIndices_[right] ].bmax_[axis]){
                --right;
            }
            if(right<=left){
                mid = left;
                break;
            }
            swap(primitiveIndices_[left], primitiveIndices_[right]);
            ++left;
            --right;
        }
#else
        PrimitivePolicy::sort(numPrimitives, &primitiveIndices_[start], bestCentroids);
        for(s32 i=start; i<end; ++i){
            if(separate < bestCentroids[primitiveIndices_[i]]){
                mid = i;
                break;
            }
        }
#endif
        if(mid == start || mid == (end-1)){
            splitMid(axis, num_l, num_r, bbox_l, bbox_r, start, numPrimitives, nodeIndex);
        } else{

            getBBox(bbox_l, start, mid);
            getBBox(bbox_r, mid, end);

            num_l = mid - start;
            num_r = numPrimitives - num_l;
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
    {
        std::ofstream file(filename, std::ios::binary);
        if(!file.is_open()){
            return;
        }

        for(s32 i=0; i<nodes_.size(); ++i){
            file << "[" << i << "] ";
            if(nodes_[i].isLeaf()){
                file << "leaf:true, face:" << nodes_[i].getPrimitiveIndex() << ":" << nodes_[i].getNumPrimitives();
            }else{
                file << "leaf:false, child:" << nodes_[i].child_;// << ", axis:" << nodes_[i].axis_;
            }

            file << ", bbox:(" << nodes_[i].bbox_.bmin_.x_ << "," << nodes_[i].bbox_.bmin_.y_ << "," << nodes_[i].bbox_.bmin_.z_ << ") - (";
            file << nodes_[i].bbox_.bmax_.x_ << "," << nodes_[i].bbox_.bmax_.y_ << "," << nodes_[i].bbox_.bmax_.z_ << ")" << std::endl;
        }
        file.close();
    }
}
#endif //INC_LRENDER_BINBVH_H__
