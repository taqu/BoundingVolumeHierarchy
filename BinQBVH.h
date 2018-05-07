#ifndef INC_ACCEL_BINQBVH_H__
#define INC_ACCEL_BINQBVH_H__
/**
@file BinQBVH.h
@author t-sakai
@date 2018/01/22 create
*/
#include "accel.h"

namespace accel
{
    template<class PrimitiveType, class PrimitivePolicy = PrimitivePolicy<PrimitiveType> >
    class BinQBVH
    {
    public:
        static constexpr f32 Epsilon = 1.0e-6f;
        static const s32 MinLeafPrimitives = 15;
        static const s32 NumBins = 32;
        static const s32 MaxBinningDepth = 11;
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

        BinQBVH();
        ~BinQBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);
        s32 getDepth() const{ return depth_;}

        void print(const char* filename);
    private:
        BinQBVH(const BinQBVH&) = delete;
        BinQBVH& operator=(const BinQBVH&) = delete;

        static const s32 MaxWorks = MaxDepth<<2;

        inline void getBBox(AABB& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 numPrimitives, const AABB& bbox);
        void split(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, f32 invArea, s32 start, s32 numPrimitives, const AABB& bbox);
        void splitMid(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 numPrimitives, const AABB& bbox);
        void splitBinned(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, f32 area, s32 start, s32 numPrimitives, const AABB& bbox);

        f32 SAH_KI_;
        f32 SAH_KT_;
        const PrimitiveType* primitives_;

        s32 depth_;
        Array<Node> nodes_;
        Array<s32> primitiveIndices_;
        Array<f32> primitiveCentroids_;
        Array<AABB> primitiveBBoxes_;
        Work works_[MaxWorks];
    };

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::Node::setJoint(s32 child, const AABB bbox[4], u8 axis[3])
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
    BinQBVH<PrimitiveType, PrimitivePolicy>::BinQBVH()
        :SAH_KI_(1.5f)
        ,SAH_KT_(1.0f)
        ,primitives_(NULL)
        ,depth_(0)
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    BinQBVH<PrimitiveType, PrimitivePolicy>::~BinQBVH()
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
    {
        f32 depth = logf(static_cast<f32>(numPrimitives) / MinLeafPrimitives) / logf(4.0f);
        s32 numNodes = static_cast<s32>(powf(2.0f, depth) + 0.5f);
        nodes_.reserve(numNodes);
        nodes_.resize(1);

        primitives_ = primitives;
        primitiveIndices_.resize(numPrimitives);
        primitiveCentroids_.resize(numPrimitives*3);
        primitiveBBoxes_.resize(numPrimitives);

        //各primitiveのcentroid, bboxを事前計算
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
    inline void BinQBVH<PrimitiveType, PrimitivePolicy>::getBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[primitiveIndices_[i]]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 numPrimitives, const AABB& bbox)
    {
        AABB childBBox[4];
        s32 primStart[4];
        s32 num[4];
        u8 axis[4];

        s32 stack = 0;
        works_[0] = Work(0, numPrimitives, 0, 1, bbox);
        while(0<=stack){
            Work work = works_[stack];
            --stack;

            depth_ = maximum(work.depth_, depth_);
            if(work.numPrimitives_<=MinLeafPrimitives || MaxDepth<=work.depth_ || MaxNodes<=nodes_.size()){
                nodes_[work.node_].setLeaf(work.start_, work.numPrimitives_);
                continue;
            }

            primStart[0] = work.start_;
            if(MaxBinningDepth<work.depth_){
                //Split top
                splitMid(axis[0], num[0], num[2], childBBox[0], childBBox[2], primStart[0], work.numPrimitives_, work.bbox_);
                primStart[2] = work.start_ + num[0];

                //Split left
                splitMid(axis[1], num[0], num[1], childBBox[0], childBBox[1], work.start_, num[0], childBBox[0]);
                primStart[1] = work.start_ + num[0];

                //Split right
                splitMid(axis[2], num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);
                primStart[3] = primStart[2] + num[2];

            } else{
                //Split top
                f32 area = work.bbox_.halfArea();
                if(area<=Epsilon){
                    splitMid(axis[0], num[0], num[2], childBBox[0], childBBox[2], primStart[0], work.numPrimitives_, work.bbox_);

                } else if(work.numPrimitives_<NumBins){
                    splitMid(axis[0], num[0], num[2], childBBox[0], childBBox[2], primStart[0], work.numPrimitives_, work.bbox_);

                } else{
                    splitBinned(axis[0], num[0], num[2], childBBox[0], childBBox[2], area, primStart[0], work.numPrimitives_, work.bbox_);
                }
                primStart[2] = work.start_ + num[0];

                //Split left
                area = childBBox[0].halfArea();
                if(area<=Epsilon){
                    splitMid(axis[1], num[0], num[1], childBBox[0], childBBox[1], primStart[0], num[0], childBBox[0]);

                } else if(num[0]<NumBins){
                    splitMid(axis[1], num[0], num[1], childBBox[0], childBBox[1], primStart[0], num[0], childBBox[0]);

                } else{
                    splitBinned(axis[1], num[0], num[1], childBBox[0], childBBox[1], area, primStart[0], num[0], childBBox[0]);
                }
                primStart[1] = work.start_ + num[0];

                //Split right
                area = childBBox[2].halfArea();
                if(area<=Epsilon){
                    splitMid(axis[2], num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);

                } else if(num[2]<NumBins){
                    splitMid(axis[2], num[2], num[3], childBBox[2], childBBox[3], primStart[2], num[2], childBBox[2]);

                } else{
                    splitBinned(axis[2], num[2], num[3], childBBox[2], childBBox[3], area, primStart[2], num[2], childBBox[2]);
                }
                primStart[3] = primStart[2] + num[2];
            }

            if(nodes_.capacity()<(nodes_.size()+4)){
                nodes_.reserve(nodes_.capacity()<<1);
            }

            s32 child = nodes_.size();
            nodes_[work.node_].setJoint(child, childBBox, axis);
            nodes_.resize(nodes_.size()+4);
            for(s32 i=0; i<4; ++i){
                works_[++stack] = Work(primStart[i], num[i], child, work.depth_+1, childBBox[i]);
                ++child;
            }
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::split(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, f32 invArea, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        s32 end = start + numPrimitives;
        s32 mid=start+(numPrimitives >> 1);

        f32 area_l, area_r;
        f32 bestCost = std::numeric_limits<f32>::max();

        //SAH, 全ての分割を試す
        axis = static_cast<u8>(bbox.maxExtentAxis());
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
                bbox_l = bl;
                bbox_r = br;
            }
        }

        num_l = mid-start;
        num_r = numPrimitives - num_l;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void BinQBVH<PrimitiveType, PrimitivePolicy>::splitMid(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, s32 start, s32 numPrimitives, const AABB& bbox)
    {
        //最大の軸を半分に分割
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
    void BinQBVH<PrimitiveType, PrimitivePolicy>::splitBinned(u8& axis, s32& num_l, s32& num_r, AABB& bbox_l, AABB& bbox_r, f32 area, s32 start, s32 numPrimitives, const AABB& bbox)
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

        Vector3 extent =bbox.extent();
        Vector3 unit = extent * (1.0f/NumBins);
        for(u8 curAxis=0; curAxis<3; ++curAxis){
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
            for(s32 i=1; i<binRight; ++i){
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

                //それぞれの分割で反転するか. 2x2x2
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
#endif //INC_ACCEL_BINQBVH_H__
