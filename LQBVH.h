#ifndef INC_ACCEL_LQBVH_H__
#define INC_ACCEL_LQBVH_H__
/**
@file LQBVH.h
@author t-sakai
@date 2018/01/22 create
*/
#include "accel.h"

namespace accel
{
    template<class PrimitiveType, class PrimitivePolicy = PrimitivePolicy<PrimitiveType> >
    class LQBVH
    {
    public:
        static constexpr f32 Epsilon = 1.0e-6f;
        static const s32 MinLeafPrimitives = 15;
        static const s32 MaxDepth = 16;
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

        union Node
        {
            static const u8 LeafFlag = (0x01U<<7);

            bool isLeaf() const
            {
                return LeafFlag == (leaf_.flags_ & LeafFlag);
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

        LQBVH();
        ~LQBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);
        s32 getDepth() const{ return depth_;}

        void print(const char* filename);
    private:
        LQBVH(const LQBVH&) = delete;
        LQBVH& operator=(const LQBVH&) = delete;

        static u32 separateBy2(u32 x);
        static u32 mortonCode3(u32 x, u32 y, u32 z);
        static bool cmpCode(const IndexCode& x0, const IndexCode& x1);

        Vector3 calcInvUnit(const AABB& bbox);
        u32 calcMortonCode3(const Vector3& x, const Vector3& invUnit, const AABB& bbox);

        void getBBox(AABB& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 numPrimitives);
        static void findSplit(s32& split, u8& axis, const IndexCode* codes, s32 first, s32 last);
        //static s32 findSplit(const IndexCode* codes, s32 first, s32 last, u32 mask);
        //static void findSplit(s32& split, u8& axis, const IndexCode* codes, s32 first, s32 last, u32 mask, s32 levels);

        const PrimitiveType* primitives_;

        s32 depth_;
        Array<Node> nodes_;
        Array<IndexCode> indexCodes_;
        Array<AABB> primitiveBBoxes_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    void LQBVH<PrimitiveType, PrimitivePolicy>::Node::setJoint(s32 child, const AABB bbox[4], u8 axis[3])
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
    LQBVH<PrimitiveType, PrimitivePolicy>::LQBVH()
        :primitives_(NULL)
        ,depth_(0)
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    LQBVH<PrimitiveType, PrimitivePolicy>::~LQBVH()
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void LQBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
    {
        s32 depth = static_cast<s32>(ceilf(logf(static_cast<f32>(numPrimitives>>4)) / logf(4.0f)));
        depth = minimum(depth, MaxDepth);
        s32 numNodes = 2;
        s32 leaves=1;
        for(s32 i=1; i<depth; ++i){
            leaves *= 4;
            numNodes += leaves;
        }
        
        nodes_.reserve(numNodes);
        nodes_.resize(1);

        primitives_ = primitives;
        indexCodes_.resize(numPrimitives);
        primitiveBBoxes_.resize(numPrimitives);

        //Calc bbox
        AABB bbox;
        bbox.setInvalid();
        for(s32 i=0; i<numPrimitives; ++i){
            primitiveBBoxes_[i] = PrimitivePolicy::getBBox(primitives_[i]);
            bbox.extend(primitiveBBoxes_[i]);
        }

        //Calc Morton codes
        indexCodes_.resize(numPrimitives);
        Vector3 invUnit = calcInvUnit(bbox);
        for(s32 i = 0; i<numPrimitives; ++i){
            Vector3 centroid = PrimitivePolicy::getCentroid(primitives_[i]);
            indexCodes_[i].index_ = i;
            indexCodes_[i].code_ = calcMortonCode3(centroid, invUnit, bbox);
        }

        depth_ = 1;
        accel::introsort(numPrimitives, &indexCodes_[0], cmpCode);
        recursiveConstruct(numPrimitives);

        primitiveBBoxes_.clear();
    }

    template<class PrimitiveType, class PrimitivePolicy>
    u32 LQBVH<PrimitiveType, PrimitivePolicy>::separateBy2(u32 x)
    {
        x = (x | (x << 8)) & 0x0000F00FU;
        x = (x | (x << 4)) & 0x000C30C3U;
        x = (x | (x << 2)) & 0x00249249U;
        return x;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    u32 LQBVH<PrimitiveType, PrimitivePolicy>::mortonCode3(u32 x, u32 y, u32 z)
    {
        return separateBy2(x) | (separateBy2(y) << 1) | (separateBy2(z) << 2);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    bool LQBVH<PrimitiveType, PrimitivePolicy>::cmpCode(const IndexCode& x0, const IndexCode& x1)
    {
        return x0.code_<x1.code_;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    Vector3 LQBVH<PrimitiveType, PrimitivePolicy>::calcInvUnit(const AABB& bbox)
    {
        Vector3 invUnit = bbox.bmax_ - bbox.bmin_;
        invUnit = invUnit * (1.0f/NumSplits);
        for(s32 i=0; i<3; ++i){
            invUnit[i] = 1.0f/invUnit[i];
        }
        return invUnit;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    u32 LQBVH<PrimitiveType, PrimitivePolicy>::calcMortonCode3(const Vector3& x, const Vector3& invUnit, const AABB& bbox)
    {
        Vector3 d = (x - bbox.bmin_);
        d *= invUnit;

        u32 v[3];
        for(s32 i=0; i<3; ++i){
            v[i] = static_cast<u32>(d[i]);
            v[i] = (v[i]<NumSplits)? v[i] : NumSplits-1;
        }
        return mortonCode3(v[0], v[1], v[2]);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void LQBVH<PrimitiveType, PrimitivePolicy>::getBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[indexCodes_[i].index_]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void LQBVH<PrimitiveType, PrimitivePolicy>::findSplit(s32& split, u8& axis, const IndexCode* codes, s32 first, s32 last)
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

#if 0
    template<class PrimitiveType, class PrimitivePolicy>
    s32 LQBVH<PrimitiveType, PrimitivePolicy>::findSplit(const IndexCode* codes, s32 first, s32 last, u32 mask)
    {
        while((first+1) < last){
            ACC_ASSERT(first<last);
            s32 mid = (first + last) >> 1;
            if((codes[first].code_ & mask) == (codes[mid].code_ & mask)){
                first = mid;
            }else{
                last = mid;
            }
        }
        return last;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void LQBVH<PrimitiveType, PrimitivePolicy>::findSplit(s32& split, u8& axis, const IndexCode* codes, s32 first, s32 last, u32 mask, s32 levels)
    {
        s32 mind = 0x7FFFFFFF;
        for(u8 i=0; i<3; ++i){
            u32 m = mask;
            for(s32 l=0; l<=levels; ++l){
                s32 s = findSplit(codes, first, last, m);
                if(first<s || s<last){
                    split = s;
                    axis = i;
                    return;
                }
                m >>= 3;
            }
            mask <<= 1;
        }
        split = (first + last)>>1;
        axis = 0;
    }
#endif

    template<class PrimitiveType, class PrimitivePolicy>
    void LQBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 numPrimitives)
    {
        Work works[MaxDepth<<2];
        AABB childBBox[4];
        s32 primStart[4];
        s32 num[4];
        u8 axis[4];

        s32 stack = 0;
        works[0] = Work(0, numPrimitives, 0, 1);
        while(0<=stack){
            Work work = works[stack];
            --stack;

            depth_ = maximum(work.depth_, depth_);
            if(work.numPrimitives_<=MinLeafPrimitives || MaxDepth<=work.depth_){
                nodes_[work.node_].setLeaf(work.start_, work.numPrimitives_);
                continue;
            }

            s32 end = work.start_ + work.numPrimitives_;
            s32 last = end-1;

            primStart[0] = work.start_;
            findSplit(primStart[2], axis[0], &indexCodes_[0], work.start_, last);
            findSplit(primStart[1], axis[1], &indexCodes_[0], work.start_, primStart[2]);
            findSplit(primStart[3], axis[2], &indexCodes_[0], primStart[2], last);

            ++primStart[1];
            ++primStart[2];
            ++primStart[3];
            num[0] = primStart[1] - work.start_;
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
            nodes_.resize(nodes_.size()+4);
            for(s32 i=0; i<4; ++i){
                works[++stack] = Work(primStart[i], num[i], child, work.depth_+1);
                ++child;
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HitRecord LQBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
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

                //‚»‚ê‚¼‚ê‚Ì•ªŠ„‚Å”½“]‚·‚é‚©. 2x2x2
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
    void LQBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
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
#endif //INC_ACCEL_LQBVH_H__
