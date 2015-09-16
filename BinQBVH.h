#ifndef INC_LRENDER_BINQBVH_H__
#define INC_LRENDER_BINQBVH_H__
/**
@file BinQBVH.h
@author t-sakai
@date 2015/09/09 create
*/
#include "lrender.h"

namespace lrender
{
    template<class PrimitiveType, class PrimitivePolicy = lrender::PrimitivePolicy<PrimitiveType> >
    class BinQBVH
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


            inline static bool isLeaf(u32 index)
            {
                return (0 != (index&LeafMask));
            }

            inline static bool isEmpty(u32 index)
            {
                return (0 == (index&LeafPrimitiveMask));
            }

            inline static u32 getPrimitiveIndex(u32 index)
            {
                return (index >> MaxNumLeafPrimitiveShift) & MaxNumPrimitives;
            }

            inline static u32 getNumPrimitives(u32 index)
            {
                return (index&LeafPrimitiveMask);
            }

            bool isJoint(s32 index) const
            {
                return (LeafMask&children_[index]) == 0;
            }

            bool isLeaf(s32 index) const
            {
                return (LeafMask&children_[index]) != 0;
            }

            bool isEmpty(s32 index) const
            {
                return (EmptyMask == children_[index]);
            }

            void setLeaf(s32 index, u32 primitiveIndex, u32 numPrimitives)
            {
                children_[index] =  LeafMask | ((primitiveIndex&MaxNumPrimitives) << MaxNumLeafPrimitiveShift) | (numPrimitives&LeafPrimitiveMask);
            }

            void setBBox(const BBox bbox[4]);

            u32 getPrimitiveIndex(s32 index) const
            {
                return (children_[index] >> MaxNumLeafPrimitiveShift) & MaxNumPrimitives;
            }

            u32 getNumPrimitives(s32 index) const
            {
                return (children_[index]&LeafPrimitiveMask);
            }

            void clear();

            __m128 bbox_[2][3];
            u32 children_[4];
            s32 axis0_;
            s32 axis1_;
            s32 axis2_;
            s32 fill_;
        };

        BinQBVH();
        ~BinQBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);
        s32 getDepth() const{ return depth_;}

        void print(const char* filename);
    private:
        BinQBVH(const BinQBVH&);
        BinQBVH& operator=(const BinQBVH&);

        inline void getBBox(BBox& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 start, s32 numPrimitives, s32 nodeIndex, const BBox& bbox, s32 depth);
        void split(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 invArea, s32 start, s32 numPrimitives, const BBox& bbox);
        void splitMid(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, s32 start, s32 numPrimitives, const BBox& bbox);
        void splitBinned(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 area, s32 start, s32 numPrimitives, const BBox& bbox);

        f32 SAH_KI_;
        f32 SAH_KT_;
        const PrimitiveType* primitives_;

        s32 depth_;
        vector_arena<Node, Align16Allocator> nodes_;
        vector_arena<s32> primitiveIndices_;
        vector_arena<f32> primitiveCentroids_;
        vector_arena<BBox> primitiveBBoxes_;
        s32* stack_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    const f32 BinQBVH<PrimitiveType, PrimitivePolicy>::Epsilon = 1.0e-6f;

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::Node::setBBox(const BBox bbox[4])
    {
        LIME_ALIGN16 f32 bb[2][3][4];

        for(s32 i=0; i<3; ++i){
            for(s32 j=0; j<4; ++j){
                bb[0][i][j] = bbox[j].bmin_[i] - Epsilon;
                bb[1][i][j] = bbox[j].bmax_[i] + Epsilon;
            }
        }

        for(s32 i=0; i<2; ++i){
            for(s32 j=0; j<3; ++j){
                _mm_store_ps((f32*)&bbox_[i][j], _mm_load_ps(bb[i][j]));
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::Node::clear()
    {
        __m128 flt_min = _mm_set1_ps(-FLT_MAX);
        __m128 flt_max = _mm_set1_ps(FLT_MAX);
        for(s32 i=0; i<3; ++i){
            _mm_store_ps((f32*)&bbox_[0][i], flt_max);
            _mm_store_ps((f32*)&bbox_[1][i], flt_min);
        }
        for(s32 i=0; i<4; ++i){
            children_[i] = EmptyMask;
        }
        axis0_ = axis1_ = axis2_ = 0;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    BinQBVH<PrimitiveType, PrimitivePolicy>::BinQBVH()
        :SAH_KI_(1.0f)
        ,SAH_KT_(1.0f)
        ,primitives_(NULL)
        ,depth_(0)
        ,stack_(NULL)
    {
        nodes_.resize(1);
        nodes_[0].clear();
    }

    template<class PrimitiveType, class PrimitivePolicy>
    BinQBVH<PrimitiveType, PrimitivePolicy>::~BinQBVH()
    {
        LIME_FREE(stack_);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
    {
        f32 depth = logf(static_cast<f32>(numPrimitives) / MinLeafPrimitives) / logf(4.0f);
        s32 numNodes = static_cast<s32>(powf(2.0f, depth) + 0.5f);
        nodes_.reserve(numNodes);
        nodes_.resize(1);
        nodes_[0].clear();

        primitives_ = primitives;
        primitiveIndices_.resize(numPrimitives);
        primitiveCentroids_.resize(numPrimitives*3);
        primitiveBBoxes_.resize(numPrimitives);

        //各primitiveのcentroid, bboxを事前計算
        f32* centroidX = &primitiveCentroids_[0];
        f32* centroidY = centroidX + numPrimitives;
        f32* centroidZ = centroidY + numPrimitives;

        BBox bbox;
        bbox.setInvalid();
        for(s32 i=0; i<numPrimitives; ++i){
            primitiveIndices_[i] = i;

            Vector3 centroid = PrimitivePolicy::getCentroid(primitives_[i]);
            centroidX[i] = centroid.x_;
            centroidY[i] = centroid.y_;
            centroidZ[i] = centroid.z_;
            primitiveBBoxes_[i] = PrimitivePolicy::getBBox(primitives_[i]);

            bbox.extend( PrimitivePolicy::getBBox(primitives_[i]) );
        }

        s32 prevDepth = depth_;
        depth_ = 1;
        if(numPrimitives<=MinLeafPrimitives){
            BBox childBBox[4];
            childBBox[0] = bbox;
            childBBox[1].setInvalid();
            childBBox[2].setInvalid();
            childBBox[3].setInvalid();

            nodes_[0].setBBox(childBBox);
            nodes_[0].setLeaf(0, 0, numPrimitives);

        }else{
            recursiveConstruct(0, numPrimitives, 0, bbox, 1);
        }

        primitiveCentroids_.clear();
        primitiveBBoxes_.clear();

        //intersectのためにスタック準備
        if(prevDepth<depth_){
            LIME_FREE(stack_);
            stack_ = (s32*)LIME_MALLOC(sizeof(s32)*depth_);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    inline void BinQBVH<PrimitiveType, PrimitivePolicy>::getBBox(BBox& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[primitiveIndices_[i]]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 start, s32 numPrimitives, s32 nodeIndex, const BBox& bbox, s32 depth)
    {
        BBox childBBox[4];
        s32 primStart[4];
        s32 num[4];
        {
            Node& node = nodes_[nodeIndex];
            primStart[0] = start;

            if(MaxBinningDepth<depth){
                splitMid(node.axis0_, num[0], num[2], childBBox[0], childBBox[2], primStart[0], numPrimitives, bbox);
                primStart[2] = start + num[0];

                splitMid(node.axis1_, num[0], num[1], childBBox[0], childBBox[1], primStart[0], num[0], childBBox[0]);
                primStart[1] = start + num[0];

                splitMid(node.axis2_, num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);
                primStart[3] = primStart[2] + num[2];

            } else{
                //top分割
                f32 area = bbox.halfArea();
                if(area<=Epsilon){
                    splitMid(node.axis0_, num[0], num[2], childBBox[0], childBBox[2], primStart[0], numPrimitives, bbox);

                } else if(numPrimitives<NumBins){
                    split(node.axis0_, num[0], num[2], childBBox[0], childBBox[2], 1.0f/area, primStart[0], numPrimitives, bbox);

                } else{
                    splitBinned(node.axis0_, num[0], num[2], childBBox[0], childBBox[2], area, primStart[0], numPrimitives, bbox);
                }
                primStart[2] = start + num[0];

                //left分割
                area = childBBox[0].halfArea();
                if(area<=Epsilon){
                    splitMid(node.axis1_, num[0], num[1], childBBox[0], childBBox[1], primStart[0], num[0], childBBox[0]);

                } else if(num[0]<NumBins){
                    split(node.axis1_, num[0], num[1], childBBox[0], childBBox[1], 1.0f/area, primStart[0], num[0], childBBox[0]);

                } else{
                    splitBinned(node.axis1_, num[0], num[1], childBBox[0], childBBox[1], area, primStart[0], num[0], childBBox[0]);
                }
                primStart[1] = start + num[0];

                //right分割
                area = childBBox[2].halfArea();
                if(area<=Epsilon){
                    splitMid(node.axis2_, num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);

                } else if(num[2]<NumBins){
                    split(node.axis2_, num[2], num[3], childBBox[2], childBBox[3], 1.0f/area, primStart[2], num[2], childBBox[2]);

                } else{
                    splitBinned(node.axis2_, num[2], num[3], childBBox[2], childBBox[3], area, primStart[2], num[2], childBBox[2]);
                }
                primStart[3] = primStart[2] + num[2];
            }

            node.setBBox(childBBox);

            for(s32 i=0; i<4; ++i){
                if(num[i]<=MinLeafPrimitives){
                    nodes_[nodeIndex].setLeaf(i, primStart[i], num[i]);
                }else{
                    s32 childIndex = nodes_.size();
                    nodes_.push_back(Node());
                    nodes_[childIndex].clear();
                    nodes_[nodeIndex].children_[i] = childIndex;
                }
            }
        }

        ++depth;
        for(s32 i=0; i<4; ++i){
            if(nodes_[nodeIndex].isJoint(i)){
                recursiveConstruct(primStart[i], num[i], nodes_[nodeIndex].children_[i], childBBox[i], depth);
            }
        }
        depth_ = maximum(depth, depth_);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::split(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 invArea, s32 start, s32 numPrimitives, const BBox& bbox)
    {
        s32 end = start + numPrimitives;
        s32 mid=start+(numPrimitives >> 1);

        f32 area_l, area_r;
        f32 bestCost = std::numeric_limits<f32>::max();
#if 1
        //SAH, 全ての分割を試す
        axis = bbox.maxExtentAxis();
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
    void BinQBVH<PrimitiveType, PrimitivePolicy>::splitMid(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, s32 start, s32 numPrimitives, const BBox& bbox)
    {
        //最大の軸を半分に分割
        axis = bbox.maxExtentAxis();

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
    void BinQBVH<PrimitiveType, PrimitivePolicy>::splitBinned(s32& axis, s32& num_l, s32& num_r, BBox& bbox_l, BBox& bbox_r, f32 area, s32 start, s32 numPrimitives, const BBox& bbox)
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

        Vector3 extent =bbox.extent();
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
                ++minBins[minIndex];
                ++maxBins[maxIndex];
            }

            Vector3 e = extent; e[curAxis] = unit[curAxis];
            f32 unitArea = e.halfArea();

            s32 binLeft = 0;
            s32 binRight = NumBins - 1;

            while(minBins[binLeft]<=0){++binLeft;}
            while(maxBins[binRight]<=0){--binRight;}

<<<<<<< HEAD
            s32 n_l = 0;
            s32 n_r = maxBins[binRight];
            for(s32 i=0; i<binRight; ++i){
=======
            s32 n_l = minBins[0];
            s32 n_r = maxBins[0];
            for(s32 i=1; i<=binRight; ++i){
>>>>>>> 5881e07bf1c0abe8ff4f64b028aa15e9a80b1243
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

        f32 separate = unit[axis] * (midBin+1) + bbox.bmin_[axis];
        s32 mid = start+(numPrimitives >> 1);

        s32 left = start;
        s32 right = end-1;

        for(;;){
            while(bestCentroids[ primitiveIndices_[left] ] < separate){
                ++left;
            }
            while(separate <= bestCentroids[ primitiveIndices_[right] ]){
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

        if(mid == start || mid == (end-1)){
            splitMid(axis, num_l, num_r, bbox_l, bbox_r, start, numPrimitives, bbox);
        } else{

            getBBox(bbox_l, start, mid);
            getBBox(bbox_r, mid, end);

            num_l = mid - start;
            num_r = numPrimitives - num_l;
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HitRecord BinQBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
    {
        __m128 origin[3];
        __m128 invDir[3];
        __m128 tminSSE;
        __m128 tmaxSSE;
        origin[0] = _mm_set1_ps(ray.origin_.x_);
        origin[1] = _mm_set1_ps(ray.origin_.y_);
        origin[2] = _mm_set1_ps(ray.origin_.z_);

        invDir[0] = _mm_set1_ps(ray.invDirection_.x_);
        invDir[1] = _mm_set1_ps(ray.invDirection_.y_);
        invDir[2] = _mm_set1_ps(ray.invDirection_.z_);

        tminSSE = _mm_set1_ps(HitEpsilon);
        tmaxSSE = _mm_set1_ps(ray.t_);

        s32 raySign[3];
        raySign[0] = (0.0f<=ray.direction_[0])? 0 : 1;
        raySign[1] = (0.0f<=ray.direction_[1])? 0 : 1;
        raySign[2] = (0.0f<=ray.direction_[2])? 0 : 1;

        HitRecord hitRecord;
        hitRecord.t_ = ray.t_;
        hitRecord.primitive_ = NULL;

        s32 stack = 0;
        u32 nodeStack[64];
        nodeStack[0] = 0;
        while(0<=stack){
            u32 index = nodeStack[stack];
            --stack;
            if(Node::isLeaf(index)){
                if(Node::isEmpty(index)){
                    continue;
                }
                u32 primIndex = Node::getPrimitiveIndex(index);
                u32 primEnd = primIndex + Node::getNumPrimitives(index);
                for(u32 i=primIndex; i<primEnd; ++i){
                    f32 t;
                    s32 idx = primitiveIndices_[i];
                    if(!primitives_[idx].testRay(t, ray)){
                        continue;
                    }
                    if(HitEpsilon < t && t < hitRecord.t_){
                        ray.t_ = t;
                        hitRecord.t_ = t;
                        hitRecord.primitive_ = &primitives_[idx];
                        tmaxSSE = _mm_set1_ps(t);
                    }
                }//for(u32 i=primIndex;

            }else{
                const Node& node = nodes_[index];
                s32 hit = qbvh::testRayAABB(tminSSE, tmaxSSE, origin, invDir, raySign, node.bbox_);

                s32 split = raySign[node.axis0_] + (raySign[node.axis1_]<<1) + (raySign[node.axis2_]<<2);

                //それぞれの分割で反転するか. 2x2x2
                static const u16 TraverseOrder[] =
                {
                    0x0123U, 0x2301U, 0x1023U, 0x3201U, 0x0132U, 0x2301U, 0x1032U, 0x3210U,
                };

                u16 order = TraverseOrder[split];

                for(s32 i=0; i<4; ++i){
                    u16 o = order&0x03U;
                    if(hit&(0x01U<<o)){
                        nodeStack[++stack] = node.children_[o];
                    }
                    order >>= 4;
                }
            }
        }//while(0<=stack){
        return hitRecord;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
    {
        std::ofstream file(filename, std::ios::binary);
        if(!file.is_open()){
            return;
        }

        for(s32 i=0; i<nodes_.size(); ++i){
            file << "[" << i << "]\n";

            const f32* bminx = (f32*)(&nodes_[i].bbox_[0][0]);
            const f32* bminy = (f32*)(&nodes_[i].bbox_[0][1]);
            const f32* bminz = (f32*)(&nodes_[i].bbox_[0][2]);

            const f32* bmaxx = (f32*)(&nodes_[i].bbox_[1][0]);
            const f32* bmaxy = (f32*)(&nodes_[i].bbox_[1][1]);
            const f32* bmaxz = (f32*)(&nodes_[i].bbox_[1][2]);

            for(s32 j=0; j<4; ++j){
                if(nodes_[i].isLeaf(j)){
                    file << "  leaf: " << nodes_[i].getPrimitiveIndex(j) << ":" << nodes_[i].getNumPrimitives(j);
                }else{
                    file << "  joint: " << nodes_[i].children_[j];
                }

                file << ", bbox:(" << bminx[j] << "," << bminy[j] << "," << bminz[j] << ") - (";
                file << bmaxx[j] << "," << bmaxy[j] << "," << bmaxz[j] << ")" << std::endl;
            }
        }
        file.close();
    }
}
#endif //INC_LRENDER_BINQBVH_H__
