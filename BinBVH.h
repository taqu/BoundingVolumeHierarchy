#ifndef INC_ACCEL_BINBVH_H__
#define INC_ACCEL_BINBVH_H__
/**
@file BinBVH.h
@author t-sakai
@date 2018/01/22 create
*/
#include "accel.h"
#include <fstream>

namespace accel
{
    template<class PrimitiveType, class PrimitivePolicy = PrimitivePolicy<PrimitiveType> >
    class BinBVH
    {
    public:
        static constexpr f32 Epsilon = 1.0e-6f;
        static const s32 MinLeafPrimitives = 15;
        static const s32 NumBins = 32;
        static const s32 MaxBinningDepth = 11;
        static const s32 MaxDepth = 24;

        struct Joint
        {
            AABB bbox_l_;
            AABB bbox_r_;
            u32 flags_;
        };

        struct Leaf
        {
            u32 padding_[10];
            s32 start_;
            s32 size_;
            u32 flags_;
        };

        union Node
        {
            static const u32 LeafFlag = ~(((u32)-1)>>1);

            Node()
            {}

            bool isLeaf() const
            {
                return LeafFlag == (leaf_.flags_ & LeafFlag);
            }

            void setLeaf(u32 start, u32 size)
            {
                leaf_.flags_ = LeafFlag;
                leaf_.start_ = start;
                leaf_.size_ = size;
            }

            void setJoint(const AABB& bbox_l, const AABB& bbox_r, s32 child)
            {
                ACC_ASSERT(0<=child);
                joint_.bbox_l_ = bbox_l;
                joint_.bbox_r_ = bbox_r;
                joint_.flags_ = child;
            }

            s32 getChildIndex() const
            {
                return joint_.flags_ & ~LeafFlag;
            }

            u32 getPrimitiveIndex() const
            {
                return leaf_.start_;
            }

            u32 getNumPrimitives() const
            {
                return leaf_.size_;
            }

            Joint joint_;
            Leaf leaf_;
        };

        struct Work
        {
            Work()
            {}

            Work(s32 start, s32 numPrimitives, s32 node, s32 depth, const AABB& bbox)
                :start_(start)
                ,numPrimitives_(numPrimitives)
                ,node_(node)
                ,depth_(depth)
                ,bbox_(bbox)
            {}

            s32 start_;
            s32 numPrimitives_;
            s32 node_;
            s32 depth_;

            AABB bbox_;
            s32 num_l_;
            s32 num_r_;
        };

        BinBVH();
        ~BinBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);

        s32 getDepth() const{ return depth_;}
        void print(const char* filename);
    private:
        BinBVH(const BinBVH&) = delete;
        BinBVH& operator=(const BinBVH&) = delete;

        inline void getBBox(AABB& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 numPrimitives, const AABB& bbox);
        void split(s32& axis, s32& num_l, s32& num_r, f32 invArea, s32 start, s32 numPrimitives, const AABB& bbox);
        void splitMid(s32& axis, s32& num_l, s32& num_r, s32 start, s32 numPrimitives, const AABB& bbox);
        void splitBinned(s32& axis, s32& num_l, s32& num_r, f32 area, s32 start, s32 numPrimitives, const AABB& bbox);

        f32 SAH_KI_;
        f32 SAH_KT_;
        const PrimitiveType* primitives_;

        s32 depth_;
        Array<Node> nodes_;
        Array<s32> primitiveIndices_;
        Array<f32> primitiveCentroids_;
        Array<AABB> primitiveBBoxes_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    BinBVH<PrimitiveType, PrimitivePolicy>::BinBVH()
        :SAH_KI_(1.5f)
        ,SAH_KT_(1.0f)
        ,primitives_(NULL)
        ,depth_(0)
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    BinBVH<PrimitiveType, PrimitivePolicy>::~BinBVH()
    {
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
        AABB bbox;
        bbox.setInvalid();

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
            bbox.extend( primitiveBBoxes_[i] );
        }

        depth_ = 1;
        recursiveConstruct(numPrimitives, bbox);

        primitiveCentroids_.clear();
        primitiveBBoxes_.clear();
    }

    template<class PrimitiveType, class PrimitivePolicy>
    inline void BinBVH<PrimitiveType, PrimitivePolicy>::getBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[primitiveIndices_[i]]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 numPrimitives, const AABB& bbox)
    {
        Work works[MaxDepth];

        s32 stack = 0;
        works[0] = Work(0, numPrimitives, 0, 1, bbox);
        s32 num_l,num_r;
        while(0<=stack){
            Work work = works[stack];
            --stack;
            depth_ = maximum(work.depth_, depth_);
            if(work.numPrimitives_ <= MinLeafPrimitives || MaxDepth<=work.depth_){
                nodes_[work.node_].setLeaf(work.start_, work.numPrimitives_);
                continue;
            }

            f32 area = work.bbox_.halfArea();
            s32 axis = 0;

#if 0
            split(axis, num_l, num_r, 1.0f/area, work.start_, work.numPrimitives_, work.bbox_);
#else
            if(MaxBinningDepth<work.depth_ || area<=Epsilon){
                splitMid(axis, num_l, num_r, work.start_, work.numPrimitives_, work.bbox_);
            } else if(work.numPrimitives_<NumBins){
                split(axis, num_l, num_r, 1.0f/area, work.start_, work.numPrimitives_, work.bbox_);
            }else{
                splitBinned(axis, num_l, num_r, area, work.start_, work.numPrimitives_, work.bbox_);
            }
#endif

            s32 childIndex = nodes_.size();
            {
                AABB bbox_l, bbox_r;
                getBBox(bbox_l, work.start_, work.start_+num_l);
                getBBox(bbox_r, work.start_+num_l, work.start_+work.numPrimitives_);
                nodes_[work.node_].setJoint(bbox_l, bbox_r, childIndex);
            }

            nodes_.resize(nodes_.size()+2);

            works[++stack] = Work(work.start_, num_l, childIndex, work.depth_+1, nodes_[work.node_].joint_.bbox_l_);
            works[++stack] = Work(work.start_+num_l, num_r, childIndex+1, work.depth_+1, nodes_[work.node_].joint_.bbox_r_);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HitRecord BinBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
    {
        HitRecord hitRecord;
        hitRecord.t_ = ray.t_;
        hitRecord.primitive_ = NULL;
        s32 top = 0;

        s32 currentNode = 0;
        s32 stack[MaxDepth<<1];
        for(;;){
            Node& node = nodes_[currentNode];

            if(node.isLeaf()){
                s32 index = node.getPrimitiveIndex();
                s32 num = node.getNumPrimitives();

                for(s32 i=0; i<num; ++i){
                    f32 t;
                    s32 idx = primitiveIndices_[index+i];
                    if(!primitives_[idx].testRay(t, ray)){
                        continue;
                    }
                    if(F32_HITEPSILON < t && t < hitRecord.t_){
                        ray.t_ = t;
                        hitRecord.t_ = t;
                        hitRecord.primitive_ = &primitives_[idx];
                    }
                }
                if(top <= 0){
                    break;
                }
                currentNode = stack[--top];

            }else{
                f32 leftMin, leftMax;
                bool hitLeft = node.joint_.bbox_l_.testRay(leftMin, leftMax, ray);
                leftMin = maximum(F32_HITEPSILON, leftMin);
                leftMax = minimum(hitRecord.t_, leftMax);

                f32 rightMin, rightMax;
                bool hitRight =  node.joint_.bbox_r_.testRay(rightMin, rightMax, ray);
                rightMin = maximum(F32_HITEPSILON, rightMin);
                rightMax = minimum(hitRecord.t_, rightMax);

                s32 childIndex = node.getChildIndex();

                if(hitLeft && !hitRight){
                    currentNode = childIndex;
                }else if(!hitLeft && hitRight){
                    currentNode = childIndex + 1;
                }else if(hitLeft && hitRight){
                    currentNode = childIndex;
                    s32 farNode = childIndex+1;
                    if(rightMin<leftMin){
                        swap(currentNode, farNode);
                    }
                    stack[top++] = farNode;
                }else{
                    if(top<=0){
                        break;
                    }
                    currentNode = stack[--top];
                }
            }//if(node.isLeaf()){
        }//for(;;){
        return hitRecord;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::split(s32& axis, s32& num_l, s32& num_r, f32 invArea, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        s32 end = start + numPrimitives;
        s32 mid=start+(numPrimitives >> 1);

        f32 area_l, area_r;
        f32 bestCost = std::numeric_limits<f32>::max();

        //SAH, 全ての分割を試す
        axis = bbox.maxExtentAxis();
        f32* bestCentroids = &primitiveCentroids_[0] + axis*primitiveIndices_.size();
        PrimitivePolicy::insertionsort(numPrimitives, &primitiveIndices_[start], bestCentroids);

        AABB bl, br;
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
            }
        }
        num_l = mid-start;
        num_r = numPrimitives - num_l;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::splitMid(s32& axis, s32& num_l, s32& num_r, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        //最大の軸を半分に分割
        axis = bbox.maxExtentAxis();

        num_l = (numPrimitives >> 1);
        num_r = numPrimitives - num_l;

        f32* centroids = &primitiveCentroids_[0] + axis * primitiveIndices_.size();
        PrimitivePolicy::sort(numPrimitives, &primitiveIndices_[start], centroids);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinBVH<PrimitiveType, PrimitivePolicy>::splitBinned(s32& axis, s32& num_l, s32& num_r, f32 area, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        ACC_ALIGN16 s32 minBins[NumBins];
        ACC_ALIGN16 s32 maxBins[NumBins];
        __m128 zero = _mm_setzero_ps();

        f32 invArea = 1.0f/area;
        axis = 0;
        s32 end = start + numPrimitives;

        f32* centroids = &primitiveCentroids_[0];
        f32* bestCentroids = centroids;
        f32 bestCost = std::numeric_limits<f32>::max();
        s32 midBin = NumBins/2;
        s32 step = static_cast<s32>(::log10f(static_cast<f32>(numPrimitives)));

        Vector3 extent = bbox.extent();
        Vector3 unit = extent * (1.0f/NumBins);
        for(s32 curAxis=0; curAxis<3; ++curAxis){
            for(s32 i=0; i<NumBins; i+=4){
                _mm_store_ps(reinterpret_cast<f32*>(&minBins[i]), zero);
                _mm_store_ps(reinterpret_cast<f32*>(&maxBins[i]), zero);
            }
            PrimitivePolicy::sort(numPrimitives, &primitiveIndices_[start], centroids);

            f32 invUnit = (absolute(unit[curAxis])<Epsilon)? 0.0f : 1.0f/unit[curAxis];
            f32 bmin = bbox.bmin_[curAxis];

            for(s32 i = start; i < end; i+=step){
                s32 index = primitiveIndices_[i];
                s32 minIndex = minimum(static_cast<s32>(invUnit * (primitiveBBoxes_[index].bmin_[curAxis] - bmin)), NumBins-1);
                s32 maxIndex = minimum(static_cast<s32>(invUnit * (primitiveBBoxes_[index].bmax_[curAxis] - bmin)), NumBins-1);
                ACC_ASSERT(0<=minIndex && minIndex<NumBins);
                ACC_ASSERT(0<=maxIndex && maxIndex<NumBins);
                ++minBins[minIndex];
                ++maxBins[maxIndex];
            }

            Vector3 e = extent; e[curAxis] = unit[curAxis];
            f32 unitArea = e.halfArea();

            s32 binLeft = 0;
            s32 binRight = NumBins - 1;

            while(minBins[binLeft]<=0){++binLeft;}
            while(maxBins[binRight]<=0){--binRight;}
            ACC_ASSERT(0<=binLeft && binLeft<NumBins);
            ACC_ASSERT(0<=binRight && binRight<NumBins);

            s32 n_l = minBins[0];
            s32 n_r = 0;
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
                    bestCentroids = centroids;
                }

                ACC_ASSERT(0<=m && m<NumBins);
                n_l += minBins[m];
                n_r -= maxBins[m];
            }

            centroids += primitiveIndices_.size();
        }//for(s32 curAxis=0;

        f32 separate = unit[axis] * (midBin+1) + bbox.bmin_[axis];
        s32 mid = start+(numPrimitives >> 1);
#if 1
        s32 left = start;
        s32 right = end-1;
        for(;;){
            while(left<end && bestCentroids[primitiveIndices_[left]]<=separate){
                ++left;
            }
            while(start<=right && separate<bestCentroids[primitiveIndices_[right]]){
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

        if(mid <= start || end<=mid){
            splitMid(axis, num_l, num_r, start, numPrimitives, bbox);
        } else{
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
#endif //INC_ACCEL_BINBVH_H__
