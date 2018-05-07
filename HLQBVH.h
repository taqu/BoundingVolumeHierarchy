#ifndef INC_ACCEL_HLQBVH_H__
#define INC_ACCEL_HLQBVH_H__
/**
@file HLQBVH.h
@author t-sakai
@date 2018/02/03 create
*/
#include "accel.h"

namespace accel
{
    template<class PrimitiveType, class PrimitivePolicy = PrimitivePolicy<PrimitiveType> >
    class HLQBVH
    {
    public:
        static constexpr f32 Epsilon = 1.0e-6f;
        static const s32 MinLeafPrimitives = 15;
        static const s32 NumBins = 16;
        static const s32 MaxDepth = 16;
        static const s32 MaxTreeletDepth = 4;
        static const s32 TreeletStart = 21;
        static const s32 NumTreelets = 64;
        static const s32 NumTreeletNodes = 85;
        static const s32 NumSplits = 1024;

        struct Joint
        {
            __m128 bbox_[2][3];
            s32 children_;
            u8 axis0_;
            u8 axis1_;
            u8 axis2_;
            u8 flags_;
        };

        struct Leaf
        {
            s32 padding0_[22];
            s32 start_;
            s32 size_;
            s32 children_;
            u8 axis0_;
            u8 axis1_;
            u8 axis2_;
            u8 flags_;
        };

        struct Treelet
        {
            s32 start_;
            s32 size_;
            AABB bbox_;
            Vector3 centroid_;
        };

        union Node
        {
            static const u8 LeafFlag = (0x01U<<7);
            static const u8 TreeletRootFlag = (0x01U<<6);

            bool isLeaf() const
            {
                return LeafFlag == (leaf_.flags_ & LeafFlag);
            }

            bool isTreeletRoot() const
            {
                return TreeletRootFlag == (leaf_.flags_ & TreeletRootFlag);
            }

            void setLeaf(s32 start, s32 size)
            {
                leaf_.flags_ = LeafFlag;
                leaf_.start_ = start;
                leaf_.size_ = size;
            }

            void setJoint(s32 child, const AABB bbox[4], u8 axis[3]);

            void setTreeletRoot(s32 start, s32 size)
            {
                leaf_.flags_ = TreeletRootFlag;
                leaf_.start_ = start;
                leaf_.size_ = size;
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

        struct IndexCode
        {
            s32 index_;
            u32 code_;
        };

        struct Work
        {
            Work()
            {}

            Work(s32 start, s32 numPrimitives, s32 node, s32 depth)
                :start_(start)
                ,numPrimitives_(numPrimitives)
                ,node_(node)
                ,depth_(depth)
            {}

            s32 start_;
            s32 numPrimitives_;
            s32 node_;
            s32 depth_;
        };

        struct TreeletWork
        {
            TreeletWork()
            {}

            TreeletWork(s32 start, s32 numPrimitives, s32 node, s32 depth, const AABB& bbox)
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
        };

        HLQBVH();
        ~HLQBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);
        s32 getDepth() const{ return depth_;}

        void print(const char* filename);
    private:
        HLQBVH(const HLQBVH&) = delete;
        HLQBVH& operator=(const HLQBVH&) = delete;

        struct SortFuncTreelet
        {
            SortFuncTreelet(u8 axis)
                :axis_(axis)
            {}

            bool operator()(const Treelet& x0, const Treelet& x1) const
            {
                f32 c0 = (x0.bbox_.bmax_[axis_] - x0.bbox_.bmin_[axis_])*0.5f;
                f32 c1 = (x1.bbox_.bmax_[axis_] - x1.bbox_.bmin_[axis_])*0.5f;
                return c0 < c1;
            }
            u8 axis_;
        };

        static u32 separateBy2(u32 x);
        static u32 mortonCode3(u32 x, u32 y, u32 z);
        static bool cmpCode(const IndexCode& x0, const IndexCode& x1);

        Vector3 calcInvUnit(const AABB& bbox);
        u32 calcMortonCode3(const Vector3& x, const Vector3& invUnit, const AABB& bbox);

        void getBBox(AABB& bbox, s32 start, s32 end);
        void getBBox(const IndexCode* indexCodes, AABB& bbox, s32 start, s32 end);
        void getTreeletBBox(AABB& bbox, s32 start, s32 end);

        void buildTreelets(const AABB& bbox);
        void recursiveConstruct(s32 start, s32 size, s32 node);
        static void findSplit(s32& split, u8& axis, const IndexCode* codes, s32 first, s32 last);

        void splitMid(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 numPrimitives, const AABB& bbox);
        void splitBinned(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, f32 area, s32 start, s32 numPrimitives, const AABB& bbox);

        struct CmpTreelet
        {
            CmpTreelet(s32 axis)
                :axis_(axis)
            {}

            bool operator()(const Treelet& x0, const Treelet& x1)
            {
                return x0.centroid_[axis_]<x1.centroid_[axis_];
            }
            s32 axis_;
        };

        const PrimitiveType* primitives_;

        f32 SAH_KI_;

        s32 depth_;
        Array<Node> nodes_;
        Array<IndexCode> indexCodes_;
        Array<AABB> primitiveBBoxes_;
        Array<Treelet> treelets_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::Node::setJoint(s32 child, const AABB bbox[4], u8 axis[3])
    {
        ACC_ALIGN16 f32 bb[2][3][4];
        joint_.flags_  = 0;
        joint_.children_ = child;
        for(s32 i=0; i<3; ++i){
            for(s32 j=0; j<4; ++j){
                bb[0][i][j] = bbox[j].bmin_[i] - Epsilon;
                bb[1][i][j] = bbox[j].bmax_[i] + Epsilon;
            }
        }

        for(s32 i=0; i<2; ++i){
            for(s32 j=0; j<3; ++j){
                _mm_store_ps((f32*)&joint_.bbox_[i][j], _mm_load_ps(bb[i][j]));
            }
        }

        joint_.axis0_ = axis[0];
        joint_.axis1_ = axis[1];
        joint_.axis2_ = axis[2];
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HLQBVH<PrimitiveType, PrimitivePolicy>::HLQBVH()
        :SAH_KI_(1.5f)
        ,primitives_(NULL)
        ,depth_(0)
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HLQBVH<PrimitiveType, PrimitivePolicy>::~HLQBVH()
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
    {
        s32 depth = static_cast<s32>(ceilf(logf(static_cast<f32>(numPrimitives>>4)) / logf(4.0f))) + MaxTreeletDepth;
        depth = minimum(depth, MaxDepth);
        s32 numNodes = 2;
        s32 leaves=1;
        for(s32 i=1; i<depth; ++i){
            leaves <<= 2;
            numNodes += leaves;
        }
        
        nodes_.reserve(numNodes);

        primitives_ = primitives;
        indexCodes_.resize(numPrimitives<<1);
        primitiveBBoxes_.resize(numPrimitives);

        //Calc bbox
        AABB bbox;
        bbox.setInvalid();
        for(s32 i=0; i<numPrimitives; ++i){
            primitiveBBoxes_[i] = PrimitivePolicy::getBBox(primitives_[i]);
            bbox.extend(primitiveBBoxes_[i]);
        }

        //Calc Morton codes
        IndexCode* indexCodes = reinterpret_cast<IndexCode*>(&indexCodes_[0]) + numPrimitives;
        Vector3 invUnit = calcInvUnit(bbox);
        for(s32 i = 0; i<numPrimitives; ++i){
            Vector3 centroid = PrimitivePolicy::getCentroid(primitives_[i]);
            indexCodes[i].index_ = i;
            indexCodes[i].code_ = calcMortonCode3(centroid, invUnit, bbox);
        }
        accel::introsort(numPrimitives, indexCodes, cmpCode);

        //Create Treelets
        static const u32 Mask = 0x3FFC0000U;
        treelets_.reserve(NumTreelets);
        treelets_.clear();
        AABB tmpBBox;
        s32 start=0;
        s32 end=1;
        for(; end<numPrimitives; ++end){
            if((indexCodes[start].code_ & Mask) != (indexCodes[end].code_ & Mask)){
                getBBox(indexCodes, tmpBBox, start, end);
                treelets_.push_back({start, end-start, tmpBBox});
                start = end;
            }
        }
        getBBox(indexCodes, tmpBBox, start, end);
        treelets_.push_back({start, end-start, tmpBBox});

        for(s32 i=0; i<treelets_.size(); ++i){
            Vector3 centroid = (treelets_[i].bbox_.bmin_ + treelets_[i].bbox_.bmax_)*0.5f;
            treelets_[i].centroid_ = centroid;
        }

        //Create SAH BVH for treelets
        buildTreelets(bbox);

        s32 count = 0;
        for(s32 i=0; i<NumTreelets; ++i){
            Node& node = nodes_[TreeletStart+i];
            ACC_ASSERT(node.isTreeletRoot());
            start = count;
            for(s32 j=0; j<node.leaf_.size_; ++j){
                const Treelet& treelet = treelets_[node.leaf_.start_+j];
                for(s32 k=0; k<treelet.size_; ++k){
                    indexCodes_[count] = indexCodes[treelet.start_+k];
                    ++count;
                }
            }
            node.leaf_.start_ = start;
            node.leaf_.size_ = count-start;
        }

        depth_ = 1;
        for(s32 i=0; i<NumTreelets; ++i){
            Node& node = nodes_[TreeletStart+i];
            ACC_ASSERT(node.isTreeletRoot());
            recursiveConstruct(node.leaf_.start_, node.leaf_.size_, TreeletStart+i);
        }

        primitiveBBoxes_.clear();
    }

    template<class PrimitiveType, class PrimitivePolicy>
    u32 HLQBVH<PrimitiveType, PrimitivePolicy>::separateBy2(u32 x)
    {
        x = (x | (x << 8)) & 0x0000F00FU;
        x = (x | (x << 4)) & 0x000C30C3U;
        x = (x | (x << 2)) & 0x00249249U;
        return x;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    u32 HLQBVH<PrimitiveType, PrimitivePolicy>::mortonCode3(u32 x, u32 y, u32 z)
    {
        return separateBy2(x) | (separateBy2(y) << 1) | (separateBy2(z) << 2);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    bool HLQBVH<PrimitiveType, PrimitivePolicy>::cmpCode(const IndexCode& x0, const IndexCode& x1)
    {
        return x0.code_<x1.code_;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    Vector3 HLQBVH<PrimitiveType, PrimitivePolicy>::calcInvUnit(const AABB& bbox)
    {
        Vector3 invUnit = bbox.bmax_ - bbox.bmin_;
        invUnit = invUnit * (1.0f/NumSplits);
        for(s32 i=0; i<3; ++i){
            invUnit[i] = 1.0f/invUnit[i];
        }
        return invUnit;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    u32 HLQBVH<PrimitiveType, PrimitivePolicy>::calcMortonCode3(const Vector3& x, const Vector3& invUnit, const AABB& bbox)
    {
        Vector3 d = (x - bbox.bmin_);
        d *= invUnit;

        u32 v[3];
        for(s32 i=0; i<3; ++i){
            v[i] = static_cast<u32>(d[i]);
            v[i] = (NumSplits<=v[i])? NumSplits-1 : v[i];
        }
        return mortonCode3(v[0], v[1], v[2]);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::getBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[indexCodes_[i].index_]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::getBBox(const IndexCode* indexCodes, AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[indexCodes[i].index_]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::getTreeletBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(treelets_[i].bbox_);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::findSplit(s32& split, u8& axis, const IndexCode* codes, s32 first, s32 last)
    {
        u32 firstCode = codes[first].code_;
        u32 lastCode = codes[last].code_;
        if(firstCode == lastCode){
            axis = 0;
            split = (first+last)>>1;
            return;
        }
        u32 commonPrefix = leadingzero(firstCode^lastCode);
        split = first;
        s32 step = last-first;
        //binary search
        do{
            step = (step+1)>>1;
            s32 newSplit = split + step;
            if(newSplit<last){
                u32 splitCode = codes[newSplit].code_;
                u32 splitPrefix = leadingzero(firstCode^splitCode);
                if(commonPrefix<splitPrefix){
                    split = newSplit;
                }
            }
        }while(1<step);
        ACC_ASSERT(first<=split && split<last);
        if(first<split){
            firstCode = codes[split-1].code_;
            lastCode = codes[split].code_;
            u32 bitIndex = 31-leadingzero(firstCode^lastCode);
            axis = bitIndex%3;
        } else{
            axis = 0;
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::buildTreelets(const AABB& bbox)
    {
        AABB childBBox[4];
        s32 primStart[4];
        s32 num[4];
        u8 axis[3];

        TreeletWork works[NumTreeletNodes];
        s32 index = 0;
        nodes_.resize(1);
        works[0] = TreeletWork(0, treelets_.size(), 0, 1, bbox);
        s32 numWorks = 1;
        while(index<NumTreeletNodes){
            TreeletWork& work = works[index];
            ++index;
            if(MaxTreeletDepth<=work.depth_){
                nodes_[work.node_].setTreeletRoot(work.start_, work.numPrimitives_);
                continue;
            }

            primStart[0] = work.start_;
            //Split top
#if 1
            f32 area = work.bbox_.halfArea();
            if(area<=Epsilon){
                splitMid(axis[0], num[0], num[2], childBBox[0], childBBox[2], primStart[0], work.numPrimitives_, work.bbox_);
            } else{
                splitBinned(axis[0], num[0], num[2], childBBox[0], childBBox[2], area, primStart[0], work.numPrimitives_, work.bbox_);
            }
#else
            splitMid(axis[0], num[0], num[2], childBBox[0], childBBox[2], primStart[0], work.numPrimitives_, work.bbox_);
#endif
            primStart[2] = work.start_ + num[0];

            //Split left
#if 1
            area = childBBox[0].halfArea();
            if(area<=Epsilon){
                splitMid(axis[1], num[0], num[1], childBBox[0], childBBox[1], primStart[0], num[0], childBBox[0]);
            } else{
                splitBinned(axis[1], num[0], num[1], childBBox[0], childBBox[1], area, primStart[0], num[0], childBBox[0]);
            }
#else
            splitMid(axis[1], num[0], num[1], childBBox[0], childBBox[1], primStart[0], num[0], childBBox[0]);
#endif
            primStart[1] = work.start_ + num[0];

            //Split right
#if 1
            area = childBBox[2].halfArea();
            if(area<=Epsilon){
                splitMid(axis[2], num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);
            } else{
                splitBinned(axis[2], num[2], num[3], childBBox[2], childBBox[3], area, primStart[2], num[2], childBBox[2]);
            }
#else
            splitMid(axis[2], num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);
#endif
            primStart[3] = primStart[2] + num[2];

            ACC_ASSERT(work.node_<nodes_.size());
            s32 child = nodes_.size();
            nodes_[work.node_].setJoint(child, childBBox, axis);
            nodes_.resize(nodes_.size()+4);
            for(s32 i=0; i<4; ++i){
                works[numWorks++] = TreeletWork(primStart[i], num[i], child, work.depth_+1, childBBox[i]);
                ++child;
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 start, s32 size, s32 node)
    {

        AABB childBBox[4];
        s32 primStart[4];
        s32 num[4];
        u8 axis[3];
        Work works[MaxDepth<<2];

        s32 stack = 0;
        works[0] = Work(start, size, node, MaxTreeletDepth);
        while(0<=stack){
            Work work = works[stack];
            --stack;

            depth_ = maximum(work.depth_, depth_);
            if(work.numPrimitives_<=MinLeafPrimitives || MaxDepth<=work.depth_){
                nodes_[work.node_].setLeaf(work.start_, work.numPrimitives_);
                continue;
            }

            start = work.start_;
            s32 end = work.start_ + work.numPrimitives_;
            s32 last = end-1;

            primStart[0] = start;
            findSplit(primStart[2], axis[0], &indexCodes_[0], start, last);
            findSplit(primStart[1], axis[1], &indexCodes_[0], start, primStart[2]);
            findSplit(primStart[3], axis[2], &indexCodes_[0], primStart[2], last);

            ++primStart[1];
            ++primStart[2];
            ++primStart[3];
            num[0] = primStart[1] - primStart[0];
            num[1] = primStart[2] - primStart[1];
            num[2] = primStart[3] - primStart[2];
            num[3] = end - primStart[3];

            getBBox(childBBox[0], primStart[0], primStart[1]);
            getBBox(childBBox[1], primStart[1], primStart[2]);
            getBBox(childBBox[2], primStart[2], primStart[3]);
            getBBox(childBBox[3], primStart[3], end);

            if(nodes_.capacity()<(nodes_.size()+4)){
                nodes_.reserve(nodes_.capacity()<<1);
            }

            s32 child = nodes_.size();
            nodes_[work.node_].setJoint(child, childBBox, axis);
            nodes_.resize(child+4);
            for(s32 i=0; i<4; ++i){
                works[++stack] = Work(primStart[i], num[i], child, work.depth_+1);
                ++child;
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::splitMid(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        //ç≈ëÂÇÃé≤Çîºï™Ç…ï™äÑ
        axis = static_cast<u8>(bbox.maxExtentAxis());

        s32 end = start + numPrimitives;
        num_l = (numPrimitives >> 1);
        num_r = numPrimitives - num_l;
        s32 mid=start+num_l;

        accel::introsort(numPrimitives, &treelets_[start], SortFuncTreelet(axis));

        getTreeletBBox(bbox_l, start, mid);
        getTreeletBBox(bbox_r, mid, end);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::splitBinned(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, f32 area, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        ACC_ALIGN16 s32 minBins[NumBins];
        ACC_ALIGN16 s32 maxBins[NumBins];

        __m128 zero = _mm_setzero_ps();

        f32 invArea = 1.0f/area;
        axis = 0;
        s32 end = start + numPrimitives;

        f32 bestCost = std::numeric_limits<f32>::max();
        s32 midBin = NumBins/2;

        Vector3 extent =bbox.extent();
        Vector3 unit = extent * (1.0f/NumBins);
        for(u8 curAxis=0; curAxis<3; ++curAxis){
            for(s32 i=0; i<NumBins; i+=4){
                _mm_store_ps(reinterpret_cast<f32*>(&minBins[i]), zero);
                _mm_store_ps(reinterpret_cast<f32*>(&maxBins[i]), zero);
            }
            accel::insertionsort(numPrimitives, &treelets_[start], SortFuncTreelet(curAxis));

            f32 invUnit = (unit[curAxis]<Epsilon)? 0.0f : 1.0f/unit[curAxis];
            f32 bmin = bbox.bmin_[curAxis];

            for(s32 i = start; i < end; ++i){
                const AABB& b = treelets_[i].bbox_;
                s32 minIndex = minimum(static_cast<s32>(invUnit * (b.bmin_[curAxis] - bmin)), NumBins-1);
                s32 maxIndex = minimum(static_cast<s32>(invUnit * (b.bmax_[curAxis] - bmin)), NumBins-1);
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
                f32 cost = SAH_KI_*invArea*(area_l*n_l + area_r*n_r);
                if(cost<bestCost){
                    midBin = m;
                    bestCost = cost;
                    axis = curAxis;
                }

                ACC_ASSERT(0<=m && m<NumBins);
                n_l += minBins[m];
                n_r -= maxBins[m];
            }
        }//for(s32 curAxis=0;

        f32 separate = unit[axis] * (midBin+1) + bbox.bmin_[axis];
        s32 mid = start+(numPrimitives >> 1);

#if 1
        s32 left = start;
        s32 right = end-1;
        for(;;){
            while(left<end && treelets_[left].centroid_[axis]<=separate){
                ++left;
            }
            while(start<=right && separate<treelets_[right].centroid_[axis]){
                --right;
            }
            if(right<=left){
                mid = left;
                break;
            }
            swap(treelets_[left], treelets_[right]);
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
            splitMid(axis, num_l, num_r, bbox_l, bbox_r, start, numPrimitives, bbox);
        } else{

            getTreeletBBox(bbox_l, start, mid);
            getTreeletBBox(bbox_r, mid, end);

            num_l = mid - start;
            num_r = numPrimitives - num_l;
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HitRecord HLQBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
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

        tminSSE = _mm_set1_ps(F32_HITEPSILON);
        tmaxSSE = _mm_set1_ps(ray.t_);

        s32 raySign[3];
        raySign[0] = (0.0f<=ray.direction_[0])? 0 : 1;
        raySign[1] = (0.0f<=ray.direction_[1])? 0 : 1;
        raySign[2] = (0.0f<=ray.direction_[2])? 0 : 1;

        HitRecord hitRecord;
        hitRecord.t_ = ray.t_;
        hitRecord.primitive_ = NULL;

        s32 stack = 0;
        s32 nodeStack[MaxDepth<<2];
        nodeStack[0] = 0;
        while(0<=stack){
            s32 index = nodeStack[stack];
            const Node& node = nodes_[index];
            ACC_ASSERT(node.leaf_.flags_ == node.joint_.flags_);
            --stack;
            if(node.isLeaf()){
                u32 primIndex = node.getPrimitiveIndex();
                u32 primEnd = primIndex + node.getNumPrimitives();
                for(u32 i=primIndex; i<primEnd; ++i){
                    f32 t;
                    s32 idx = indexCodes_[i].index_;
                    if(!primitives_[idx].testRay(t, ray)){
                        continue;
                    }
                    if(F32_HITEPSILON < t && t < hitRecord.t_){
                        ray.t_ = t;
                        hitRecord.t_ = t;
                        hitRecord.primitive_ = &primitives_[idx];
                        tmaxSSE = _mm_set1_ps(t);
                    }
                }//for(u32 i=primIndex;

            }else{
                s32 hit = qbvh::testRayAABB(tminSSE, tmaxSSE, origin, invDir, raySign, node.joint_.bbox_);
                s32 split = raySign[node.joint_.axis0_] + (raySign[node.joint_.axis1_]<<1) + (raySign[node.joint_.axis2_]<<2);

                //ÇªÇÍÇºÇÍÇÃï™äÑÇ≈îΩì]Ç∑ÇÈÇ©. 2x2x2
                static const u16 TraverseOrder[] =
                {
                    0x0123U, 0x2301U, 0x1023U, 0x3201U, 0x0132U, 0x2301U, 0x1032U, 0x3210U,
                };
                u16 order = TraverseOrder[split];
                s32 children = node.joint_.children_;
                for(s32 i=0; i<4; ++i){
                    u16 o = order&0x03U;
                    if(hit&(0x01U<<o)){
                        nodeStack[++stack] = children + o;
                    }
                    order >>=4;
                }
            }
        }//while(0<=stack){
        return hitRecord;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void HLQBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
    {
        std::ofstream file(filename, std::ios::binary);
        if(!file.is_open()){
            return;
        }

        for(s32 i=0; i<nodes_.size(); ++i){
            file << "[" << i << "]\n";

            for(s32 j=0; j<4; ++j){
                if(nodes_[i].isLeaf()){
                    file << "  leaf: " << nodes_[i].getPrimitiveIndex() << ":" << nodes_[i].getNumPrimitives();

                    const f32* bminx = (f32*)(&nodes_[i].joint_.bbox_[0][0]);
                    const f32* bminy = (f32*)(&nodes_[i].joint_.bbox_[0][1]);
                    const f32* bminz = (f32*)(&nodes_[i].joint_.bbox_[0][2]);

                    const f32* bmaxx = (f32*)(&nodes_[i].joint_.bbox_[1][0]);
                    const f32* bmaxy = (f32*)(&nodes_[i].joint_.bbox_[1][1]);
                    const f32* bmaxz = (f32*)(&nodes_[i].joint_.bbox_[1][2]);

                    file << ", bbox:(" << bminx[j] << "," << bminy[j] << "," << bminz[j] << ") - (";
                    file << bmaxx[j] << "," << bmaxy[j] << "," << bmaxz[j] << ")" << std::endl;

                }else{
                    file << "  joint: " << std::endl;
                }
            }
        }
        file.close();
    }
}
#endif //INC_ACCEL_HLQBVH_H__
