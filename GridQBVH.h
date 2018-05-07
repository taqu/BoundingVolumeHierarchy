#ifndef INC_ACCEL_GRIDQBVH_H_
#define INC_ACCEL_GRIDQBVH_H_
/**
@file GRIDQBVH.h
@author t-sakai
@date 2018/02/08 create
*/
#include "accel.h"
#ifdef _DEBUG
#include <stdio.h>

#endif

namespace accel
{
    template<class PrimitiveType, class PrimitivePolicy = PrimitivePolicy<PrimitiveType> >
    class GRIDQBVH
    {
    public:
        static constexpr f32 Epsilon = 1.0e-6f;
        static const s32 MinLeafPrimitives = 8;
        static const s32 MaxGridDepth = 9;
        static const s32 MaxLastLevel = 3;
        static const s32 MaxDepth = MaxGridDepth+MaxLastLevel;

        static const s32 TopGridShift = 5;
        static const s32 SecondGridShift = 5;
        static const s32 TotalShift = TopGridShift+SecondGridShift;
        static const s32 CellSplits = 1<<TotalShift;
        static const s32 TopGridResolution = (1<<TopGridShift);
        static const s32 SecondGridResolution = (1<<SecondGridShift);

        static const s32 TopCodeShift = TopGridShift*3;
        static const s32 TopCodeMask = (1<<TopCodeShift)-1;
        static const s32 SecondCodeShift = SecondGridShift*3;
        static const s32 SecondCodeMask = (1<<SecondCodeShift)-1;

        static const s32 IntMax = std::numeric_limits<s32>::max();
        static const s32 IntMin = std::numeric_limits<s32>::min();

        static const u8 TopFlag = (0x01U<<5);
        static const u8 SecondFlag = (0x01U<<6);
        static const u8 LeafFlag = (0x01U<<7);

        struct Range
        {
            s32 start_;
            s32 size_;
        };

        struct TriReference
        {
            static const u32 TopMask = (0x01U<<TopGridShift)-1;
            static const u32 SecondMask = (0x01U<<SecondGridShift)-1;

            void getTop(s32& x, s32& y, s32& z) const
            {
                accel::rmortonCode3(x, y, z, cellId_);
                x = (x>>SecondGridShift) & TopMask;
                y = (y>>SecondGridShift) & TopMask;
                z = (z>>SecondGridShift) & TopMask;
            }

            void getSecond(s32& x, s32& y, s32& z) const
            {
                accel::rmortonCode3(x, y, z, cellId_);
                x &= SecondMask;
                y &= SecondMask;
                z &= SecondMask;
            }

            s32 id_;
            s32 cellId_;
        };

        struct Grid
        {
            s32 start_;
            s32 size_;
            u8 top_[4];
        };

        struct IAABB
        {
            s32 min_[3];
            s32 max_[3];

            s32 maxExtentAxis() const
            {
                s32 size[3] = {max_[0]-min_[0], max_[1]-min_[1], max_[2]-min_[2]};
                s32 axis = (size[0] < size[1])
                    ? (size[1] < size[2])? 2 : 1
                    : (size[0] < size[2])? 2 : 0;
                return axis;
            }
        };

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
            s32 padding0_[21];
            s32 start_;
            s32 size_;
            s32 level_;
            s32 children_;
            u8 axis0_;
            u8 axis1_;
            u8 axis2_;
            u8 flags_;
        };

        union Node
        {
            bool isTop() const
            {
                return TopFlag == (leaf_.flags_ & TopFlag); 
            }

            bool isSecond() const
            {
                return SecondFlag == (leaf_.flags_ & SecondFlag);
            }

            bool isLeaf() const
            {
                return LeafFlag == (leaf_.flags_ & LeafFlag);
            }

            void setTop(s32 start, s32 size)
            {
                leaf_.flags_ = TopFlag;
                leaf_.start_ = start;
                leaf_.size_ = size;
            }

            void setSecond(s32 start, s32 size, s32 level)
            {
                leaf_.flags_ = SecondFlag;
                leaf_.start_ = start;
                leaf_.size_ = size;
                leaf_.level_ = level;
            }

            void setLeaf(s32 start, s32 size)
            {
                leaf_.flags_ = LeafFlag;
                leaf_.start_ = start;
                leaf_.size_ = size;
            }

            void setJoint(s32 child, const AABB bbox[4], u8 axis[3]);

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

        union BBox
        {
            BBox()
            {}

            BBox(const AABB& bbox)
                :fbbox_(bbox)
            {}
            BBox(const IAABB& ibbox)
                :ibbox_(ibbox)
            {}

            IAABB ibbox_;
            AABB fbbox_;
        };

        struct Work
        {
            Work()
            {}

            Work(s32 start, s32 size, s32 node, s32 depth, const AABB& bbox)
                :start_(start)
                ,size_(size)
                ,node_(node)
                ,depth_(depth)
                ,bbox_(bbox)
            {}

            Work(s32 start, s32 size, s32 node, s32 depth, const IAABB& bbox)
                :start_(start)
                ,size_(size)
                ,node_(node)
                ,depth_(depth)
                ,bbox_(bbox)
            {}

            s32 start_;
            s32 size_;
            s32 node_;
            s32 depth_;
            BBox bbox_;
        };

        GRIDQBVH();
        ~GRIDQBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);
        s32 getDepth() const{ return depth_;}

        void print(const char* filename);
    private:
        GRIDQBVH(const GRIDQBVH&) = delete;
        GRIDQBVH& operator=(const GRIDQBVH&) = delete;

        static constexpr s32 calcTotalNodes(s32 levels)
        {
            s32 total = 1;
            while(0<--levels){
                total <<= 2;
            }
            return total;
        }

        static const s32 MaxWorks = MaxDepth<<2;
        static const s32 MaxNodes = calcTotalNodes(MaxDepth);

        struct SortFunc
        {
            SortFunc(u8 axis)
                :axis_(axis)
            {}

            bool operator()(const Grid& grid0, const Grid& grid1) const
            {
                return grid0.top_[axis_] < grid1.top_[axis_];
            }

            u8 axis_;
        };

        static void radixSort(s32 size, TriReference* tries);

        static void buildTopGrid(Array<Grid>& grid, s32 size, const TriReference* tries);
        static void buildSecondGrid(Array<Grid>& grid, Node node, const TriReference* tries);
        void buildTopHierarchy();
        void buildSecondHierarchy(s32 root);
        void buildLastHierarchy(s32 root, s32 level);

        void getBBox(IAABB& bbox, s32 start, s32 end);
        void getBBox(AABB& bbox, s32 start, s32 end);
        void getBBoxFromReferences(AABB& bbox, s32 start, s32 end);

        void splitSAH(u8& axis, s32& num_l, s32& num_r, IAABB& bbox_l, IAABB& bbox_r, s32 start, s32 size, const IAABB& bbox);
        void split(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 size, const AABB& bbox);

        Range getRange(const Grid* begin, const Grid* end);

        static void sort(u8 axis, s32 num, Grid* grids);
        static s32 divide(u8 axis, s32 shift, s32 m, s32 size, s32 bmin, Grid* grids);
        static s32 divide(u8 axis, f32 mid, s32 start, s32 size, const PrimitiveType* primitives, TriReference* tries);

        const PrimitiveType* primitives_;

        s32 depth_;
        Array<Node> nodes_;
        Array<TriReference> triReferences_;
        Array<TriReference> workReferences_;
        Array<Grid> grids_;
        Work works_[MaxWorks];
    };

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::Node::setJoint(s32 child, const AABB bbox[4], u8 axis[3])
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
    GRIDQBVH<PrimitiveType, PrimitivePolicy>::GRIDQBVH()
        :primitives_(NULL)
        ,depth_(0)
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    GRIDQBVH<PrimitiveType, PrimitivePolicy>::~GRIDQBVH()
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
    {
        nodes_.reserve(MaxNodes);
        triReferences_.resize(numPrimitives);
        workReferences_.clear();
        workReferences_.reserve(numPrimitives);
        primitives_ = primitives;

        s32 gridSize = maximum(TopGridResolution, SecondGridResolution) + 1;
        grids_.reserve(gridSize*gridSize*gridSize);

        //Calc bbox
        AABB bbox;
        bbox.setInvalid();
        for(s32 i=0; i<numPrimitives; ++i){
            AABB b = PrimitivePolicy::getBBox(primitives_[i]);
            bbox.extend(b);
        }

        //Calc Cell IDs
        Vector3 size = bbox.bmax_ - bbox.bmin_;
        Vector3 invUnit(CellSplits/size.x_, CellSplits/size.y_, CellSplits/size.z_);
        for(s32 i=0; i<numPrimitives; ++i){
            Vector3 centroid = PrimitivePolicy::getCentroid(primitives_[i]) - bbox.bmin_;
            s32 x = minimum(static_cast<s32>(centroid.x_ * invUnit.x_), CellSplits-1);
            s32 y = minimum(static_cast<s32>(centroid.y_ * invUnit.y_), CellSplits-1);
            s32 z = minimum(static_cast<s32>(centroid.z_ * invUnit.z_), CellSplits-1);
            s32 cellId = accel::mortonCode3(x, y, z);
            triReferences_[i] = {i, cellId};
        }

        //Sort by Cell IDs
        radixSort(numPrimitives, &triReferences_[0]);

        //for(s32 i=1; i<numPrimitives; ++i){
        //    ACC_ASSERT(triReferences_[i-1].cellId_<=triReferences_[i].cellId_);
        //}

        depth_ = 1;

        //Build top
        grids_.clear();
        buildTopGrid(grids_, numPrimitives, &triReferences_[0]);
        ACC_ASSERT(grids_.size()<=(gridSize*gridSize*gridSize));
        nodes_.clear();
        buildTopHierarchy();

        workReferences_.swap(triReferences_);
        workReferences_.clear();
        //Build second
        s32 secondStart = nodes_.size();
        for(s32 i=0; i<secondStart; ++i){
            if(!nodes_[i].isTop()){
                continue;
            }
            grids_.clear();
            buildSecondGrid(grids_, nodes_[i], &triReferences_[0]);
            buildSecondHierarchy(i);
        }
        workReferences_.swap(triReferences_);
        workReferences_.clear();
        grids_.clear();

        s32 secondEnd = nodes_.size();
        //Build last
        for(s32 i=secondStart; i<secondEnd; ++i){
            if(!nodes_[i].isSecond()){
                continue;
            }
            buildLastHierarchy(i, nodes_[i].leaf_.level_);
        }
    }


    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::radixSort(s32 size, TriReference* tries)
    {
        static const s32 NumPasses = 6;
        static const s32 NumBits = 30;
        static const s32 BitsPerPass = NumBits/NumPasses;
        static const s32 NumBuckets = 1<<BitsPerPass;
        static const s32 Mask = NumBuckets-1;
        s32 bucketCount[NumBuckets];
        s32 indices[NumBuckets];

        TriReference* tmp = reinterpret_cast<TriReference*>(ACC_MALLOC(sizeof(TriReference)*size));
        TriReference* in = tries;
        TriReference* out = tmp;
        for(s32 pass=0; pass<NumPasses; ++pass){
            memset(bucketCount, 0, sizeof(s32)*NumBuckets);
            s32 shift = pass*BitsPerPass;
            for(s32 i=0; i<size; ++i){
                s32 bucket = (in[i].cellId_>>shift) & Mask;
                ++bucketCount[bucket];
            }
            indices[0] = 0;
            for(s32 i=1; i<NumBuckets; ++i){
                indices[i] = indices[i-1] + bucketCount[i-1];
            }
            for(s32 i=0; i<size; ++i){
                s32 bucket = (in[i].cellId_>>shift) & Mask;
                out[indices[bucket]++] = in[i];
            }
            swap(in, out);

        }//for(s32 i=0;
        ACC_FREE(tmp);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::buildTopGrid(Array<Grid>& grid, s32 size, const TriReference* tries)
    {
        if(size<=0){
            return;
        }

        s32 x,y,z, prev;
        prev = (tries[0].cellId_>>TopCodeShift) & TopCodeMask;
        s32 start = 0;
        tries[start].getTop(x,y,z);

        for(s32 i=1; i<size; ++i){
#if _DEBUG
            s32 tx,ty,tz;
            tries[i].getTop(tx,ty,tz);
            ACC_ASSERT(tz<TopGridResolution);
            ACC_ASSERT(ty<TopGridResolution);
            ACC_ASSERT(tz<TopGridResolution);
#endif
            s32 index = (tries[i].cellId_>>TopCodeShift) & TopCodeMask;
            if(prev != index){
                ACC_ASSERT(grid.size()<grid.capacity());
                Grid g;
                g.start_ = start;
                g.size_ = i-start;
                g.top_[0] = static_cast<u8>(x);
                g.top_[1] = static_cast<u8>(y);
                g.top_[2] = static_cast<u8>(z);

                grid.push_back(g);
                prev = index;
                start = i;
                tries[start].getTop(x,y,z);
            }
        }

        Grid g;
        g.start_ = start;
        g.size_ = size-start;
        g.top_[0] = static_cast<u8>(x);
        g.top_[1] = static_cast<u8>(y);
        g.top_[2] = static_cast<u8>(z);
        grid.push_back(g);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::buildSecondGrid(Array<Grid>& grid, Node node, const TriReference* tries)
    {
        ACC_ASSERT(0<node.leaf_.size_);

        s32 start = node.leaf_.start_;
        s32 end = start + node.leaf_.size_;
        s32 x,y,z, prev;
        tries[start].getSecond(x,y,z);
        prev = tries[start].cellId_ & SecondCodeMask;

        for(s32 i=start+1; i<end; ++i){
#if _DEBUG
            s32 tx,ty,tz;
            tries[i].getSecond(tx,ty,tz);
            ACC_ASSERT(tz<SecondGridResolution);
            ACC_ASSERT(ty<SecondGridResolution);
            ACC_ASSERT(tz<SecondGridResolution);
#endif
            s32 index = tries[i].cellId_ & SecondCodeMask;
            if(prev != index){
                ACC_ASSERT(grid.size()<grid.capacity());
                Grid g;
                g.start_ = start;
                g.size_ = i-start;
                g.top_[0] = static_cast<u8>(x);
                g.top_[1] = static_cast<u8>(y);
                g.top_[2] = static_cast<u8>(z);

                grid.push_back(g);
                prev = index;
                start = i;
                tries[start].getSecond(x,y,z);
            }
        }

        Grid g;
        g.start_ = start;
        g.size_ = end-start;
        g.top_[0] = static_cast<u8>(x);
        g.top_[1] = static_cast<u8>(y);
        g.top_[2] = static_cast<u8>(z);
        grid.push_back(g);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::buildTopHierarchy()
    {
        IAABB childBBox[4];
        s32 start[4];
        s32 num[4];
        u8 axis[4];

        s32 stack = 0;
        {
            IAABB bbox;
            getBBox(bbox, 0, grids_.size());
            works_[0] = Work(0, grids_.size(), 0, 1, bbox);
            nodes_.push_back(Node());
        }
        while(0<=stack){
            Work work = works_[stack];
            --stack;

            if(work.size_<=1 || MaxGridDepth<=work.depth_){
                depth_ = maximum(work.depth_, depth_);
                const Grid* grid = &grids_[work.start_];
                Range range = getRange(grid, grid+work.size_);
                nodes_[work.node_].setTop(range.start_, range.size_);
                continue;
            }

            start[0] = work.start_;

            splitSAH(axis[0], num[0], num[2], childBBox[0], childBBox[2], start[0], work.size_, work.bbox_.ibbox_);
            start[2] = work.start_ + num[0];

            //Split left
            splitSAH(axis[1], num[0], num[1], childBBox[0], childBBox[1], work.start_, num[0], childBBox[0]);
            start[1] = work.start_ + num[0];

            //Split right
            splitSAH(axis[2], num[2], num[3], childBBox[2], childBBox[3], start[2], num[2], childBBox[2]);
            start[3] = start[2] + num[2];

            if(nodes_.capacity()<(nodes_.size()+4)){
                nodes_.reserve(nodes_.capacity()<<1);
            }

            s32 child = nodes_.size();
            {
                AABB bboxes[4];
                getBBox(bboxes[0], start[0], start[0]+num[0]);
                getBBox(bboxes[1], start[1], start[1]+num[1]);
                getBBox(bboxes[2], start[2], start[2]+num[2]);
                getBBox(bboxes[3], start[3], start[3]+num[3]);

                nodes_[work.node_].setJoint(child, bboxes, axis);
            }
            nodes_.resize(nodes_.size()+4);
            for(s32 i=0; i<4; ++i, ++child){
                if(num[i]<=0){
                    nodes_[child].setLeaf(0, 0);
                    continue;
                }

                works_[++stack] = Work(start[i], num[i], child, work.depth_+1, childBBox[i]);
                ACC_ASSERT(0<=stack && stack<MaxWorks);
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::buildSecondHierarchy(s32 root)
    {
        IAABB childBBox[4];
        s32 start[4];
        s32 num[4];
        u8 axis[4];

        s32 stack = 0;
        {
            IAABB bbox;
            getBBox(bbox, 0, grids_.size());
            works_[0] = Work(0, grids_.size(), root, TopGridShift+1, bbox);
        }
        while(0<=stack){
            Work work = works_[stack];
            --stack;

            if(work.size_<=1 || MaxGridDepth<=work.depth_){
                depth_ = maximum(work.depth_, depth_);
                const Grid* grid = &grids_[work.start_];
                Range range = getRange(grid, grid+work.size_);
                if(range.size_<=MinLeafPrimitives){
                    nodes_[work.node_].setLeaf(range.start_, range.size_);
                }else{
                    nodes_[work.node_].setSecond(range.start_, range.size_, work.depth_);
                }
                continue;
            }

            start[0] = work.start_;

            splitSAH(axis[0], num[0], num[2], childBBox[0], childBBox[2], start[0], work.size_, work.bbox_.ibbox_);
            start[2] = work.start_ + num[0];

            //Split left
            splitSAH(axis[1], num[0], num[1], childBBox[0], childBBox[1], work.start_, num[0], childBBox[0]);
            start[1] = work.start_ + num[0];

            //Split right
            splitSAH(axis[2], num[2], num[3], childBBox[2], childBBox[3], start[2], num[2], childBBox[2]);
            start[3] = start[2] + num[2];

            if(nodes_.capacity()<(nodes_.size()+4)){
                nodes_.reserve(nodes_.capacity()<<1);
            }

            s32 child = nodes_.size();
            {
                AABB bboxes[4];
                getBBox(bboxes[0], start[0], start[0]+num[0]);
                getBBox(bboxes[1], start[1], start[1]+num[1]);
                getBBox(bboxes[2], start[2], start[2]+num[2]);
                getBBox(bboxes[3], start[3], start[3]+num[3]);

                nodes_[work.node_].setJoint(child, bboxes, axis);
            }
            nodes_.resize(nodes_.size()+4);
            for(s32 i=0; i<4; ++i, ++child){
                if(num[i]<=0){
                    nodes_[child].setLeaf(0, 0);
                    continue;
                }

                works_[++stack] = Work(start[i], num[i], child, work.depth_+1, childBBox[i]);
                ACC_ASSERT(0<=stack && stack<MaxWorks);
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::buildLastHierarchy(s32 root, s32 level)
    {
        AABB childBBox[4];
        s32 start[4];
        s32 num[4];
        u8 axis[4];

        s32 stack = 0;
        {
            Node& node = nodes_[root];
            AABB bbox;
            getBBoxFromReferences(bbox, node.leaf_.start_, node.leaf_.start_+node.leaf_.size_);
            works_[0] = Work(node.leaf_.start_, node.leaf_.size_, root, level, bbox);
        }

        while(0<=stack){
            Work work = works_[stack];
            --stack;

            if(work.size_<=MinLeafPrimitives || MaxDepth<=work.depth_){
                depth_ = maximum(work.depth_, depth_);
                Node& node = nodes_[work.node_];
                node.setLeaf(work.start_, work.size_);
                continue;
            }

            start[0] = work.start_;

            split(axis[0], num[0], num[2], childBBox[0], childBBox[2], start[0], work.size_, work.bbox_.fbbox_);
            start[2] = work.start_ + num[0];

            //Split left
            split(axis[1], num[0], num[1], childBBox[0], childBBox[1], work.start_, num[0], childBBox[0]);
            start[1] = work.start_ + num[0];

            //Split right
            split(axis[2], num[2], num[3], childBBox[2], childBBox[3], start[2], num[2], childBBox[2]);
            start[3] = start[2] + num[2];

            if(nodes_.capacity()<(nodes_.size()+4)){
                nodes_.reserve(nodes_.capacity()<<1);
            }

            s32 child = nodes_.size();
            nodes_[work.node_].setJoint(child, childBBox, axis);
            nodes_.resize(nodes_.size()+4);
            for(s32 i=0; i<4; ++i, ++child){
                if(num[i]<=0){
                    nodes_[child].setLeaf(0, 0);
                    continue;
                }

                works_[++stack] = Work(start[i], num[i], child, work.depth_+1, childBBox[i]);
                ACC_ASSERT(0<=stack && stack<MaxWorks);
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::getBBox(IAABB& bbox, s32 start, s32 end)
    {
        if(end<=start){
            bbox = {0,0,0, 0,0,0};
            return;
        }

        bbox = {IntMax, IntMax, IntMax, IntMin, IntMin, IntMin};
        for(s32 i=start; i<end; ++i){
            const Grid& grid = grids_[i];
            ACC_ASSERT(0<grid.size_);
            bbox.max_[0] = maximum(static_cast<s32>(grid.top_[0]), bbox.max_[0]);
            bbox.max_[1] = maximum(static_cast<s32>(grid.top_[1]), bbox.max_[1]);
            bbox.max_[2] = maximum(static_cast<s32>(grid.top_[2]), bbox.max_[2]);

            bbox.min_[0] = minimum(static_cast<s32>(grid.top_[0]), bbox.min_[0]);
            bbox.min_[1] = minimum(static_cast<s32>(grid.top_[1]), bbox.min_[1]);
            bbox.min_[2] = minimum(static_cast<s32>(grid.top_[2]), bbox.min_[2]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::getBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            const Grid& grid = grids_[i];
            ACC_ASSERT(0<grid.size_);
            s32 gridEnd = grid.start_ + grid.size_;
            for(s32 j=grid.start_; j<gridEnd; ++j){
                s32 refId = triReferences_[j].id_;
                AABB b = PrimitivePolicy::getBBox(primitives_[refId]);
                bbox.extend(b);
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::getBBoxFromReferences(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            s32 refId = triReferences_[i].id_;
            AABB b = PrimitivePolicy::getBBox(primitives_[refId]);
            bbox.extend(b);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::splitSAH(u8& axis, s32& num_l, s32& num_r, IAABB& bbox_l, IAABB& bbox_r, s32 start, s32 size, const IAABB& bbox)
    {
        if(size<=0){
            axis = 0;
            num_l = 0;
            num_r = 0;
            bbox_l = {0,0,0, 0,0,0};
            bbox_r = {0,0,0, 0,0,0};
            return;
        }
        //static const s32 NumBins = 8;
        static const s32 NumBins = 32;
        ACC_ALIGN16 s32 numTri[NumBins];
        ACC_ALIGN16 s32 numTriL[NumBins];
        ACC_ALIGN16 s32 numTriR[NumBins];
        ACC_ALIGN16 s32 areaL[NumBins];
        ACC_ALIGN16 s32 areaR[NumBins];

        __m128 zero = _mm_setzero_ps();

        s32 mid = 0;
        s32 end = start + size;
        s32 bestCost = std::numeric_limits<s32>::max();
        s32 shift = 0;
        axis = 0;
        //Calculate SAH cost for each axis
        for(u8 currentAxis=0; currentAxis<3; ++currentAxis){
            s32 extent = bbox.max_[currentAxis] - bbox.min_[currentAxis] + 1;
            ACC_ASSERT(0<extent);

            //128, 64 32 16 8
            //static const s32 MaxShift = 4;
            static const s32 MaxShift = 2;
            s32 s;
            for(s=0; s<=MaxShift; ++s){
                if((extent>>s)<=NumBins){
                    break;
                }
            }
            //for(s=MaxShift; 0<s; --s){
            //    if((extent >> s)<=NumBins){
            //        break;
            //    }
            //}

            for(s32 i=0; i<NumBins; i+=4){
                _mm_store_ps(reinterpret_cast<f32*>(&numTri[i]), zero);
                _mm_store_ps(reinterpret_cast<f32*>(&numTriL[i]), zero);
                _mm_store_ps(reinterpret_cast<f32*>(&numTriR[i]), zero);
                _mm_store_ps(reinterpret_cast<f32*>(&areaL[i]), zero);
                _mm_store_ps(reinterpret_cast<f32*>(&areaR[i]), zero);
            }

            s32 total = 0;
            for(s32 i=start; i<end; ++i){
                const Grid& grid = grids_[i];
                s32 index = (static_cast<s32>(grid.top_[currentAxis]) - bbox.min_[currentAxis]) >> s;
                ACC_ASSERT(index<NumBins);
                numTri[index] += grid.size_;
                total += grid.size_;
            }


            numTriL[0] = numTri[0];
            numTriR[0] = total;

            s32 left = 0<numTri[0]? 0 : NumBins-1;
            s32 right = 0;
            for(s32 i=1; i<NumBins; ++i){
                s32 prev = i-1;
                numTriL[i] = numTri[i] + numTriL[prev];
                numTriR[i] = numTriR[prev] - numTriL[prev];
                if(0<numTri[i]){
                    left = minimum(left, i);
                    right = maximum(right, i);
                }
            }

            for(s32 i=0; i<NumBins; ++i){
                areaL[i] = maximum(i-left+1, 0);
                areaR[i] = maximum(right-i+1, 0);
            }

            for(s32 i=1; i<NumBins; ++i){
                s32 cost = (areaL[i-1] * numTriL[i-1] + areaR[i] * numTriR[i])/extent;
                if(cost<=bestCost){
                    bestCost = cost;
                    axis = currentAxis;
                    shift = s;
                    mid = i;
                }
            }
        }

        ACC_ASSERT(0<=mid && mid<NumBins);
        mid = divide(axis, shift, mid, size, bbox.min_[axis], &grids_[start]);
        ACC_ASSERT(0<=mid && mid<size);
        num_l = mid;
        num_r = size - num_l;

        mid += start;
        getBBox(bbox_l, start, mid);
        getBBox(bbox_r, mid, end);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::split(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 size, const AABB& bbox)
    {
        if(size<=0){
            axis = 0;
            num_l = 0;
            num_r = 0;
            bbox_l.setInvalid();
            bbox_r.setInvalid();
            return;
        }
        Vector3 extent = bbox.extent();
        axis = 0;
        f32 max = extent[0];
        for(u8 i=1; i<3; ++i){
            if(max<extent[i]){
                axis = i;
                max = extent[i];
            }
        }

        f32 mid = (bbox.bmax_[axis] + bbox.bmin_[axis])*0.5f;
        s32 m = divide(axis, mid, start, size, primitives_, &triReferences_[0]);
        ACC_ASSERT(start<=m);
        num_l = m-start;
        num_r = size-num_l;
        ACC_ASSERT(0<=num_l);
        ACC_ASSERT(0<=num_r);

        getBBoxFromReferences(bbox_l, start, start+num_l);
        getBBoxFromReferences(bbox_r, start+num_l, start+size);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HitRecord GRIDQBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
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
        u32 nodeStack[MaxDepth<<2];
        nodeStack[0] = 0;
        while(0<=stack){
            u32 index = nodeStack[stack];
            const Node& node = nodes_[index];
            ACC_ASSERT(node.leaf_.flags_ == node.joint_.flags_);
            --stack;
            if(node.isLeaf()){
                s32 start = node.getPrimitiveIndex();
                s32 end = start + node.getNumPrimitives();
                for(s32 i=start; i<end; ++i){
                    f32 t;
                    s32 refId = triReferences_[i].id_;
                    if(!primitives_[refId].testRay(t, ray)){
                        continue;
                    }
                    if(F32_HITEPSILON < t && t < hitRecord.t_){
                        ray.t_ = t;
                        hitRecord.t_ = t;
                        hitRecord.primitive_ = &primitives_[refId];
                        tmaxSSE = _mm_set1_ps(t);
                    }
                }//for(u32 i=primIndex;

            }else{
                s32 hit = qbvh::testRayAABB(tminSSE, tmaxSSE, origin, invDir, raySign, node.joint_.bbox_);
                s32 split = raySign[node.joint_.axis0_] + (raySign[node.joint_.axis1_]<<1) + (raySign[node.joint_.axis2_]<<2);

                //2x2x2 patterns of traversal
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
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::sort(u8 axis, s32 num, Grid* grids)
    {
        accel::introsort(num, grids, SortFunc(axis));
    }

    template<class PrimitiveType, class PrimitivePolicy>
    typename GRIDQBVH<PrimitiveType, PrimitivePolicy>::Range GRIDQBVH<PrimitiveType, PrimitivePolicy>::getRange(const Grid* begin, const Grid* end)
    {
        s32 size = 0;
        for(const Grid* grid = begin; grid != end; ++grid){
            size += grid->size_;
        }
        s32 start = workReferences_.size();
        workReferences_.resize(workReferences_.size() + size);
        s32 count = start;
        for(const Grid* grid = begin; grid != end; ++grid){
            for(s32 i=0; i<grid->size_; ++i){
                workReferences_[count++] = triReferences_[grid->start_+i];
            }
        }
        return {start, size};
    }

    template<class PrimitiveType, class PrimitivePolicy>
    s32 GRIDQBVH<PrimitiveType, PrimitivePolicy>::divide(u8 axis, s32 shift, s32 m, s32 size, s32 bmin, Grid* grids)
    {
        if(size<=0){
            return 0;
        }
        s32 left = 0;
        s32 right = size-1;
        for(;;){
            for(;;){
                s32 l = (static_cast<s32>(grids[left].top_[axis]) - bmin) >> shift;
                if(m<=l){
                    break;
                }
                ++left;
            }
            for(;;){
                s32 r = (static_cast<s32>(grids[right].top_[axis]) - bmin) >> shift;
                if(r<=m){
                    break;
                }
                --right;
            }
            if(right<=left){
                break;
            }
            swap(grids[left], grids[right]);
            ++left;
            --right;
        }
        return right;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    s32 GRIDQBVH<PrimitiveType, PrimitivePolicy>::divide(u8 axis, f32 mid, s32 start, s32 size, const PrimitiveType* primitives, TriReference* tries)
    {
        ACC_ASSERT(0<size);
        s32 left = start;
        s32 right = start+size-1;
        for(;;){
            while(left<right){
                f32 l = PrimitivePolicy::getCentroid(primitives[tries[left].id_], axis);
                if(mid<l){
                    break;
                }
                ++left;
            }
            while(left<right){
                f32 r = PrimitivePolicy::getCentroid(primitives[tries[right].id_], axis);
                if(r<=mid){
                    break;
                }
                --right;
            }
            if(right<=left){
                break;
            }
            swap(tries[left], tries[right]);
            ++left;
            --right;
        }
        return right;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void GRIDQBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
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
#endif //INC_ACCEL_GRIDQBVH_H_
