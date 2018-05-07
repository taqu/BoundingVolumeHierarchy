#ifndef INC_ACCEL_SORT_H__
#define INC_ACCEL_SORT_H__
/**
@file Sort.h
@author t-sakai
@date 2018/01/22 create
*/
#include "accel.h"

namespace accel
{
    template<class T>
    struct SortCompFuncType
    {
        typedef bool(*SortCmp)(const T& lhs, const T& rhs);
    };

    /**
    @return true if lhs<rhs, false if lhs>=rhs
    */
    //template<class T>
    //typedef bool(*SortCompFunc)(const T& lhs, const T& rhs);

    template<class T>
    bool less(const T& lhs, const T& rhs)
    {
        return (lhs<rhs);
    }

    /**
    @return a bucket number for radix sort
    */
    template<class T>
    s32 calcBucket(const T& x)
    {
        return x;
    }

    //------------------------------------------------
    //---
    //--- insertionsort
    //---
    //------------------------------------------------
    /**

    Uはbool operator(const T& a, const T& b) const{ return a<b;}が必要
    */
    template<class T, class U>
    void insertionsort(s32 n, T* v, U func)
    {
        for(s32 i=1; i<n; ++i){
            for(s32 j=i; 0<j; --j){
                s32 k = j-1;
                if(func(v[j], v[k])){
                    swap(v[j], v[k]);
                }else{
                    break;
                }
            }
        }
    }

    template<class T>
    void insertionsort(s32 n, T* v)
    {
        insertionsort(n, v, less<T>);
    }

    //------------------------------------------------
    //---
    //--- heapsort
    //---
    //------------------------------------------------
    /**

    Uはbool operator(const T& a, const T& b) const{ return a<b;}が必要
    */
    template<class T, class U>
    void heapsort(s32 n, T* v, U func)
    {
        ACC_ASSERT(0<=n);

        --v;
        s32 i, j;
        T x;
        for(s32 k=n>>1; k>=1; --k){
            i=k;
            x = v[k];
            while((j=i<<1)<=n){
                if(j<n && func(v[j], v[j+1])){
                    ++j;
                }

                if(!func(x, v[j])){
                    break;
                }
                v[i] = v[j];
                i = j;
            }
            v[i] = x;
        }

        while(n>1){
            x = v[n];
            v[n] = v[1];
            --n;
            i = 1;
            while((j=i<<1)<=n){
                if(j<n && func(v[j], v[j+1])){
                    ++j;
                }

                if(!func(x, v[j])){
                    break;
                }
                v[i] = v[j];
                i = j;
            }
            v[i] = x;
        }
    }

    template<class T>
    void heapsort(s32 n, T* v)
    {
        heapsort(n, v, less<T>);
    }

    //------------------------------------------------
    //---
    //--- quicksort
    //---
    //------------------------------------------------
    /**

    U: bool operator(const T& a, const T& b) const{ return a<b;}
    */
    template<class T, class U>
    void quicksort(s32 n, T* v, U func)
    {
        static const s32 SwitchN = 47;
        if(n<SwitchN){
            insertionsort(n, v, func);
            return;
        }

        s32 i0 = 0;
        s32 i1 = n-1;

        T pivot = v[ (i0+i1)>>1 ];

        for(;;){
            while(func(v[i0], pivot)){
                ++i0;
            }

            while(func(pivot, v[i1])){
                --i1;
            }

            if(i1<=i0){
                break;
            }
            lcore::swap(v[i0], v[i1]);
            ++i0;
            --i1;
        }

        if(1<i0){
            quicksort(i0, v, func);
        }

        ++i1;
        n = n-i1;
        if(1<n){
            quicksort(n, v+i1, func);
        }
    }

    template<class T>
    void quicksort(s32 n, T* v)
    {
        quicksort(n, v, less<T>);
    }

    //------------------------------------------------
    //---
    //--- introsort
    //---
    //------------------------------------------------
    /**

    Uはbool operator(const T& a, const T& b) const{ return a<b;}が必要
    */
    template<class T, class U>
    void introsort(s32 n, T* v, s32 depth, U func)
    {
        static const s32 SwitchN = 47;
        if(n<SwitchN){
            insertionsort(n, v, func);
            return;
        }
        if(depth<=0){
            heapsort(n, v, func);
            return;
        }

        s32 i0 = 0;
        s32 i1 = n-1;

        T pivot = v[ (i0+i1)>>1 ];

        for(;;){
            while(func(v[i0], pivot)){
                ++i0;
            }

            while(func(pivot, v[i1])){
                --i1;
            }

            if(i1<=i0){
                break;
            }
            swap(v[i0], v[i1]);
            ++i0;
            --i1;
        }

        --depth;
        if(1<i0){
            introsort(i0, v, depth, func);
        }

        ++i1;
        n = n-i1;
        if(1<n){
            introsort(n, v+i1, depth, func);
        }
    }

    template<class T, class U>
    void introsort(s32 n, T* v, U func)
    {
        s32 depth = 0;
        s32 t = n;
        while(1<t){
            ++depth;
            t >>= 1;
        }
        introsort(n, v, depth, func);
    }

    template<class T>
    void introsort(s32 n, T* v)
    {
        introsort(n, v, less<T>);
    }

    //------------------------------------------------
    //---
    //--- radixsort
    //---
    //------------------------------------------------
    template<class T, class U, s32 NumBuckets>
    void radixsort(s32 size, T* dst, const T* src, U calcBucket)
    {
        static_assert(0<NumBuckets, "NumBuckets should be greater than 0.");
        ACC_ASSERT(0<=size);
        ACC_ASSERT(NULL != dst);
        ACC_ASSERT(NULL != src);

        s32 bucketCount[NumBuckets] = {0};
        for(s32 i=0; i<size; ++i){
            s32 bucket = calcBucket(src[i]);
            ACC_ASSERT(0<=bucket && bucket<NumBuckets);
            ++bucketCount[bucket];
        }
        s32 indices[NumBuckets];
        indices[0] = 0;
        for(s32 i=1; i<NumBuckets; ++i){
            indices[i] = indices[i-1] + bucketCount[i-1];
        }
        for(s32 i=0; i<size; ++i){
            s32 bucket = calcBucket(src[i]);
            dst[indices[bucket]++] = src[i];
        }
    }
}
#endif //INC_ACCEL_SORT_H__
