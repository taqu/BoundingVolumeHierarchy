#ifndef INC_ACCEL_MEDBVH_H__
#define INC_ACCEL_MEDBVH_H__
/**
@file MedBVH.h
@author t-sakai
@date 2018/01/22 create
*/
#include "accel.h"
#include <fstream>

namespace accel
{
    template<class PrimitiveType, class PrimitivePolicy = PrimitivePolicy<PrimitiveType> >
    class MedBVH
    {
    public:
        static constexpr f32 Epsilon = 1.0e-6f;
        static const s32 MinLeafPrimitives = 15;
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
        };

        MedBVH();
        ~MedBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);

        s32 getDepth() const{ return depth_;}
        void print(const char* filename);
    private:
        MedBVH(const MedBVH&) = delete;
        MedBVH& operator=(const MedBVH&) = delete;

        inline void getBBox(AABB& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 numPrimitives, const AABB& bbox);

        const PrimitiveType* primitives_;

        s32 depth_;
        Array<Node> nodes_;
        Array<s32> primitiveIndices_;
        Array<f32> primitiveCentroids_;
        Array<AABB> primitiveBBoxes_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    MedBVH<PrimitiveType, PrimitivePolicy>::MedBVH()
        :primitives_(NULL)
        ,depth_(0)
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    MedBVH<PrimitiveType, PrimitivePolicy>::~MedBVH()
    {
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void MedBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
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

        //Šeprimitive‚Ìcentroid, bbox‚ðŽ–‘OŒvŽZ
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
    inline void MedBVH<PrimitiveType, PrimitivePolicy>::getBBox(AABB& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[primitiveIndices_[i]]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void MedBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 numPrimitives, const AABB& bbox)
    {
        Work works[MaxDepth];

        s32 stack = 0;
        works[0] = Work(0, numPrimitives, 0, 1, bbox);
        s32 num_l,num_r;
        AABB bbox_l, bbox_r;
        while(0<=stack){
            Work work = works[stack];
            --stack;
            depth_ = maximum(work.depth_, depth_);
            if(work.numPrimitives_ <= MinLeafPrimitives || MaxDepth<=work.depth_){
                nodes_[work.node_].setLeaf(work.start_, work.numPrimitives_);
                continue;
            }

            s32 end = work.start_ + work.numPrimitives_;

            s32 mid = work.start_ + (work.numPrimitives_ >> 1);
            s32 axis = 0;

            num_l = (work.numPrimitives_ >> 1);
            num_r = work.numPrimitives_ - num_l;

            axis = work.bbox_.maxExtentAxis();
            f32* centroids = &primitiveCentroids_[0] + axis*primitiveIndices_.size();
            PrimitivePolicy::sort(work.numPrimitives_, &primitiveIndices_[work.start_], centroids);

            getBBox(bbox_l, work.start_, mid);
            getBBox(bbox_r, mid, end);

            s32 childIndex = nodes_.size();
            nodes_[work.node_].setJoint(bbox_l, bbox_r, childIndex);

            nodes_.resize(nodes_.size()+2);
            works[++stack] = Work(work.start_, num_l, childIndex, work.depth_+1, bbox_l);
            works[++stack] = Work(work.start_+num_l, num_r, childIndex+1, work.depth_+1, bbox_r);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    HitRecord MedBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
    {
        HitRecord hitRecord;
        hitRecord.t_ = ray.t_;
        hitRecord.primitive_ = NULL;
        s32 top = 0;

        s32 currentNode = 0;
        s32 stack[MaxDepth];
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
    void MedBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
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
#endif //INC_ACCEL_MEDBVH_H__
