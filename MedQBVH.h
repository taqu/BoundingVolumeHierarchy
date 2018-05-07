#ifndef INC_ACCEL_MEDQBVH_H__
#define INC_ACCEL_MEDQBVH_H__
/**
@file MedQBVH.h
@author t-sakai
@date 2018/01/22 create
*/
#include "accel.h"

namespace accel
{
    template<class PrimitiveType, class PrimitivePolicy = PrimitivePolicy<PrimitiveType> >
    class MedQBVH
    {
    public:
        static constexpr f32 Epsilon = 1.0e-6f;
        static const s32 MinLeafPrimitives = 15;
        static const s32 MaxDepth = 24;
        static const s32 MaxNodes = 0xFFFFFF-4;

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

            Node()
            {}

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

        };

        MedQBVH();
        ~MedQBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);
        s32 getDepth() const{ return depth_;}

        void print(const char* filename);
    private:
        MedQBVH(const MedQBVH&) = delete;
        MedQBVH& operator=(const MedQBVH&) = delete;

        inline void getBBox(AABB& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 numPrimitives, const AABB& bbox);
        void splitMid(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 numPrimitives, const AABB& bbox);

        const PrimitiveType* primitives_;

        s32 depth_;
        Array<Node> nodes_;
        Array<s32> primitiveIndices_;
        Array<f32> primitiveCentroids_;
        Array<AABB> primitiveBBoxes_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    void MedQBVH<PrimitiveType, PrimitivePolicy>::Node::setJoint(s32 child, const AABB bbox[4], u8 axis[3])
    {
        ACC_ALIGN16 f32 bb[2][3][4];
        joint_.flags_ = 0;
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
    MedQBVH<PrimitiveType, PrimitivePolicy>::MedQBVH()
        :primitives_(NULL)
        ,depth_(0)
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    MedQBVH<PrimitiveType, PrimitivePolicy>::~MedQBVH()
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void MedQBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
    {
        f32 depth = logf(static_cast<f32>(numPrimitives) / MinLeafPrimitives) / logf(4.0f);
        s32 numNodes = static_cast<s32>(powf(2.0f, depth) + 0.5f);
        nodes_.reserve(numNodes);
        nodes_.resize(1);

        primitives_ = primitives;
        primitiveIndices_.resize(numPrimitives);
        primitiveCentroids_.resize(numPrimitives*3);
        primitiveBBoxes_.resize(numPrimitives);

        //äeprimitiveÇÃcentroid, bboxÇéñëOåvéZ
        f32* centroidX = &primitiveCentroids_[0];
        f32* centroidY = centroidX + numPrimitives;
        f32* centroidZ = centroidY + numPrimitives;

        AABB bbox;
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

        depth_ = 1;
        recursiveConstruct(numPrimitives, bbox);

        primitiveCentroids_.clear();
        primitiveBBoxes_.clear();
    }

    template<class PrimitiveType, class PrimitivePolicy>
    inline void MedQBVH<PrimitiveType, PrimitivePolicy>::getBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[primitiveIndices_[i]]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void MedQBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 numPrimitives, const AABB& bbox)
    {
        Work works[MaxDepth<<2];
        AABB childBBox[4];
        s32 primStart[4];
        s32 num[4];
        u8 axis[4];

        s32 stack = 0;
        works[0] = Work(0, numPrimitives, 0, 1, bbox);
        while(0<=stack){
            Work work = works[stack];
            --stack;

            depth_ = maximum(work.depth_, depth_);

            if(work.numPrimitives_ <= MinLeafPrimitives || MaxDepth<=work.depth_ || MaxNodes<=nodes_.size()){
                nodes_[work.node_].setLeaf(work.start_, work.numPrimitives_);
                continue;
            }

            primStart[0] = work.start_;

            splitMid(axis[0], num[0], num[2], childBBox[0], childBBox[2], primStart[0], work.numPrimitives_, work.bbox_);
            primStart[2] = work.start_ + num[0];

            //Split left
            splitMid(axis[1], num[0], num[1], childBBox[0], childBBox[1], work.start_, num[0], childBBox[0]);
            primStart[1] = work.start_ + num[0];

            //Split right
            splitMid(axis[2], num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);
            primStart[3] = primStart[2] + num[2];

            if(nodes_.capacity()<(nodes_.size()+4)){
                nodes_.reserve(nodes_.capacity()<<1);
            }

            s32 child = nodes_.size();
            nodes_[work.node_].setJoint(child, childBBox, axis);
            nodes_.resize(nodes_.size()+4);
            for(s32 i=0; i<4; ++i){
                works[++stack] = Work(primStart[i], num[i], child, work.depth_+1, childBBox[i]);
                ++child;
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void MedQBVH<PrimitiveType, PrimitivePolicy>::splitMid(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        //ç≈ëÂÇÃé≤Çîºï™Ç…ï™äÑ
        axis = static_cast<u8>(bbox.maxExtentAxis());

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
    HitRecord MedQBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
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
                    s32 idx = primitiveIndices_[i];
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
    void MedQBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
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
#endif //INC_ACCEL_MEDQBVH_H__
