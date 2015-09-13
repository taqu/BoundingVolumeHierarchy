#ifndef INC_LRENDER_MIDBVH_H__
#define INC_LRENDER_MIDBVH_H__
/**
@file MidBVH.h
@author t-sakai
@date 2015/09/09 create
*/
#include "lrender.h"
#include <fstream>

namespace lrender
{
    template<class PrimitiveType, class PrimitivePolicy = lrender::PrimitivePolicy<PrimitiveType> >
    class MidBVH
    {
    public:
        static const f32 Epsilon;
        static const s32 MinLeafPrimitives = 15;

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

        MidBVH();
        ~MidBVH();

        void build(s32 numPrimitives, const PrimitiveType* primitives);
        HitRecord intersect(Ray& ray);

        s32 getDepth() const{ return depth_;}
        void print(const char* filename);
    private:
        MidBVH(const MidBVH&);
        MidBVH& operator=(const MidBVH&);

        inline void getBBox(BBox& bbox, s32 start, s32 end);

        void recursiveConstruct(s32 start, s32 numPrimitives, s32 nodeIndex, s32 depth);

        const PrimitiveType* primitives_;

        s32 depth_;
        vector_arena<Node> nodes_;
        vector_arena<s32> primitiveIndices_;
        vector_arena<f32> primitiveCentroids_;
        vector_arena<BBox> primitiveBBoxes_;
        s32* stack_;
    };

    template<class PrimitiveType, class PrimitivePolicy>
    const f32 MidBVH<PrimitiveType, PrimitivePolicy>::Epsilon = 1.0e-6f;

    template<class PrimitiveType, class PrimitivePolicy>
    MidBVH<PrimitiveType, PrimitivePolicy>::MidBVH()
        :primitives_(NULL)
        ,depth_(0)
        ,stack_(NULL)
    {
        nodes_.resize(1);
        nodes_[0].bbox_.zero();
        nodes_[0].child_ = Node::EmptyMask;
    }

    template<class PrimitiveType, class PrimitivePolicy>
    MidBVH<PrimitiveType, PrimitivePolicy>::~MidBVH()
    {
        LIME_FREE(stack_);
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void MidBVH<PrimitiveType, PrimitivePolicy>::build(s32 numPrimitives, const PrimitiveType* primitives)
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
    inline void MidBVH<PrimitiveType, PrimitivePolicy>::getBBox(BBox& bbox, s32 start, s32 end)
    {
        bbox.setInvalid();
        for(s32 i=start; i<end; ++i){
            bbox.extend(primitiveBBoxes_[primitiveIndices_[i]]);
        }
    }

    template<class PrimitiveType, class PrimitivePolicy>
    void MidBVH<PrimitiveType, PrimitivePolicy>::recursiveConstruct(s32 start, s32 numPrimitives, s32 nodeIndex, s32 depth)
    {
        if (numPrimitives <= MinLeafPrimitives){
            nodes_[nodeIndex].setLeaf(start, numPrimitives);
            depth_ = maximum(depth, depth_);
            return;
        }

        s32 end = start + numPrimitives;

        s32 num_l, num_r, mid=start+(numPrimitives >> 1);
        BBox bbox_l, bbox_r;
        s32 axis = 0;

        num_l = (numPrimitives >> 1);
        num_r = numPrimitives - num_l;

        axis = nodes_[nodeIndex].bbox_.maxExtentAxis();
        f32* centroids = &primitiveCentroids_[0] + axis*primitiveIndices_.size();
        PrimitivePolicy::sort(numPrimitives, &primitiveIndices_[start], centroids);

        getBBox(bbox_l, start, mid);
        getBBox(bbox_r, mid, end);


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
    HitRecord MidBVH<PrimitiveType, PrimitivePolicy>::intersect(Ray& ray)
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
    void MidBVH<PrimitiveType, PrimitivePolicy>::print(const char* filename)
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
#endif //INC_LRENDER_MIDBVH_H__
