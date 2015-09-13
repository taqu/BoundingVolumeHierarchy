#ifndef INC_LRENDER_SORT_H__
#define INC_LRENDER_SORT_H__
/**
@file Sort.h
@author t-sakai
@date 2011/12/31 create
*/
#include "lrender.h"

namespace lrender
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
        LASSERT(0<=n);

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

    Uはbool operator(const T& a, const T& b) const{ return a<b;}が必要
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
    //--- Timsort
    //---
    //------------------------------------------------
    template<class T, class Comp=SortCompFuncType<T> >
    class Timsort
    {
    public:
        static const s32 MAX_MERGE_PENDING = 85;
        static const s32 MIN_GALLOP = 7;
        static const s32 MERGESTATE_TEMP_SIZE = 256;

        typedef T value_type;
        typedef value_type* pointer_type;
        typedef const value_type* const_pointer_type;

        typedef value_type& reference_type;
        typedef const value_type& const_reference_type;
        typedef Comp comp_type;

        Timsort(comp_type compFunc)
            :compFunc_(compFunc)
        {}

        ~Timsort()
        {
            merge_free();
        }

        void sort(pointer_type keys, s32 size);
    private:
        Timsort(const Timsort&);
        Timsort& operator=(const Timsort&);

        inline static s32 minimum(s32 v0, s32 v1)
        {
            return (v0<v1)? v0 : v1;
        }

        inline static s32 maximum(s32 v0, s32 v1)
        {
            return (v0<v1)? v1 : v0;
        }

        struct sortslice
        {
            pointer_type keys_;

            static inline void copy(sortslice& s0, s32 i, sortslice& s1, s32 j)
            {
                s0.keys_[i] = s1.keys_[j];
            }

            static inline void copy_incr(sortslice& s0, sortslice& s1)
            {
                *s0.keys_++ = *s1.keys_++;
            }

            static inline void copy_decr(sortslice& s0, sortslice& s1)
            {
                *s0.keys_-- = *s1.keys_--;
            }

            static inline void memcpy(sortslice& s0, s32 i, sortslice& s1, s32 j, s32 n)
            {
                ::memcpy(&s0.keys_[i], &s1.keys_[j], sizeof(value_type)*n);
            }

            static inline void memmove(sortslice& s0, s32 i, sortslice& s1, s32 j, s32 n)
            {
                ::memmove(&s0.keys_[i], &s1.keys_[j], sizeof(value_type)*n);
            }

            static inline void advance(sortslice& slice, s32 n)
            {
                slice.keys_ += n;
            }

            static void reverse(pointer_type lo, pointer_type hi);
            static inline void reverse(sortslice& slice, s32 n)
            {
                reverse(slice.keys_, &slice.keys_[n]);
            }
        };

        struct s_slice
        {
            sortslice base_;
            s32 length_;
        };

        /**
        @brief 挿入位置を2分探索する挿入ソート
        */
        void binarysort(sortslice lo, pointer_type hi, pointer_type start);

        /**
        @brief 先頭から順序が連続する部分配列の長さを返す.
        @param lo ... 先頭
        @param hi ... 終端
        @param descending ... 出力. 降順か
        */
        s32 count_run(pointer_type lo, pointer_type hi, s32& descending);

        /**
        @return p ... a[p-1]< key <=a[p]
        @param key ...
        @param a ... ソート済み配列
        @param n ... 0<n
        @param hint ... 0<= hint <n
        */
        s32 gallop_left(const_reference_type key, pointer_type a, s32 n, s32 hint);

        /**
        @return p ... a[p-1]<= key <a[p]
        @param key ...
        @param a ... ソート済み配列
        @param n ... 0<n
        @param hint ... 0<= hint <n
        */
        s32 gallop_right(const_reference_type key, pointer_type a, s32 n, s32 hint);

        void merge_init(s32 size);
        void merge_free();
        s32 merge_getmem(s32 need);
        s32 merge_compute_minrun(s32 n);

#define MERGE_GETMEM(NEED) ((NEED)<=allocated_? 0 : merge_getmem(NEED))

        /**
        @brief ssaとssbをマージ
        ssb.keys_ == ssa.keys_+na, 0<na, 0<nb, na<=nb
        */
        void merge_lo(sortslice ssa, s32 na, sortslice ssb, s32 nb);

        /**
        @brief ssaとssbをマージ
        ssb.keys_ == ssa.keys_+na, 0<na, 0<nb, nb<=na
        */
        void merge_hi(sortslice ssa, s32 na, sortslice ssb, s32 nb);

        /**
        @brief

        pendingスタックが以下の不等式を満たすまで処理.
        1. len[-3] > len[-2] + len[-1]
        2. len[-2] > len[-1]
        */
        void merge_collapse();
        void merge_force_collapse();

        /**
        @brief スタックのn, n+1番目の部分列をマージ
        */
        void merge_at(s32 n);

        comp_type compFunc_;

        s32 min_gallop_;
        sortslice a_;
        s32 allocated_;
        s32 numPendings_;
        s_slice pending_[MAX_MERGE_PENDING];
        value_type temparray_[MERGESTATE_TEMP_SIZE];
    };

    template<class T, class Comp>
    void Timsort<T,Comp>::sortslice::reverse(pointer_type lo, pointer_type hi)
    {
        LASSERT(NULL != lo);
        LASSERT(NULL != hi);
        --hi;
        while(lo<hi){
            value_type t = *lo;
            *lo = *hi;
            *hi = t;
            ++lo;
            --hi;
        }
    }

    // 挿入位置を2分探索する挿入ソート
    template<class T, class Comp>
    void Timsort<T, Comp>::binarysort(sortslice lo, pointer_type hi, pointer_type start)
    {
        LASSERT(lo.keys_<=start && start<=hi);
        // [lo, start)はソート済み
        if(lo.keys_ == start){
            ++start;
        }

        pointer_type l;
        pointer_type p;
        pointer_type r;

        value_type pivot;

        for(; start<hi; ++start){
            l = lo.keys_;
            r = start;
            pivot = *r;

            //lo.keys_[l] == pivotとなるlを探索
            do{
                p = l + ((r-l)>>1);
                if(compFunc_(pivot, *p)){
                    r = p;
                } else{
                    l = p+1;
                }
            } while(l<r);
            LASSERT(l==r);

            //配列シフト
            for(p=start; l<p; --p){
                *p = *(p-1);
            }
            *l = pivot;
        }
    }

    // 先頭から順序が連続する部分配列の長さを返す.
    template<class T, class Comp>
    s32 Timsort<T, Comp>::count_run(pointer_type lo, pointer_type hi, s32& descending)
    {
        LASSERT(lo<hi);
        descending = 0;
        ++lo;
        if(lo == hi){
            return 1;
        }
        s32 n = 2;
        if(compFunc_(*lo, *(lo-1))){
            descending = 1;
            for(lo=lo+1; lo<hi; ++lo, ++n){
                if(false == compFunc_(*lo, *(lo-1))){
                    break;
                }
            }
        } else{
            for(lo=lo+1; lo<hi; ++lo, ++n){
                if(compFunc_(*lo, *(lo-1))){
                    break;
                }
            }
        }
        return n;
    }

    // p ... a[p-1]< key <=a[p]
    template<class T, class Comp>
    s32 Timsort<T, Comp>::gallop_left(const_reference_type key, pointer_type a, s32 n, s32 hint)
    {
        LASSERT(NULL != a);
        LASSERT(0<n);
        LASSERT(0<=hint && hint<n);

        a += hint;
        s32 lastoffset = 0;
        s32 offset = 1;
        if(compFunc_(*a, key)){
            //a[hint] < key, gallop right
            const s32 maxoffset = n-hint;
            while(offset<maxoffset){
                if(compFunc_(a[offset], key)){
                    lastoffset = offset;
                    offset = (offset<<1)+1;
                    if(offset<=0){//int overflow
                        offset = maxoffset;
                    }
                } else{
                    break;
                }
            }// while(offset<maxoffset)
            // a[hint + lastoffset] < key <= a[hint + offset]
            offset = minimum(offset, maxoffset);
            lastoffset += hint;
            offset += hint;

        } else{
            //key <= a[hint], gallop left
            const s32 maxoffset = hint+1;
            while(offset<maxoffset){
                if(compFunc_(a[-offset], key)){
                    break;
                } else{
                    lastoffset = offset;
                    offset = (offset<<1)+1;
                    if(offset<=0){//int overflow
                        offset = maxoffset;
                    }
                }
            }
            // a[hint - offset] < key <= a[hint - lastoffset]
            offset = minimum(offset, maxoffset);

            //lastoffsetとoffset入れ換え
            s32 t = lastoffset;
            lastoffset = hint - offset;
            offset = hint - t;
        }
        a -= hint;
        LASSERT(-1<=lastoffset && lastoffset<offset);
        LASSERT(offset<=n);

        //a[lastoffset] < key <= a[offset]
        //binary search
        ++lastoffset;
        while(lastoffset<offset){
            s32 m = lastoffset + ((offset-lastoffset)>>1);
            if(compFunc_(a[m], key)){
                lastoffset = m + 1;
            } else{
                offset = m;
            }
        }
        LASSERT(lastoffset == offset);
        return offset;
    }

    // p ... a[p-1]<= key <a[p]
    template<class T, class Comp>
    s32 Timsort<T, Comp>::gallop_right(const_reference_type key, pointer_type a, s32 n, s32 hint)
    {
        LASSERT(NULL != a);
        LASSERT(0<n);
        LASSERT(0<=hint && hint<n);

        a += hint;
        s32 lastoffset = 0;
        s32 offset = 1;
        if(compFunc_(key, *a)){
            //key < a[hint], gallop left
            const s32 maxoffset = hint + 1;
            while(offset<maxoffset){
                if(compFunc_(key, a[-offset])){
                    lastoffset = offset;
                    offset = (offset<<1)+1;
                    if(offset<=0){//int overflow
                        offset = maxoffset;
                    }
                } else{
                    break;
                }
            }// while(offset<maxoffset)
            // a[hint - offset] <= key < a[hint - lastoffset]
            offset = minimum(offset, maxoffset);

            //lastoffsetとoffset入れ換え
            s32 t = lastoffset;
            lastoffset = hint - offset;
            offset = hint - t;

        } else{
            //key <= a[hint], gallop left
            const s32 maxoffset = n-hint;
            while(offset<maxoffset){
                if(compFunc_(key, a[offset])){
                    break;
                } else{
                    lastoffset = offset;
                    offset = (offset<<1)+1;
                    if(offset<=0){//int overflow
                        offset = maxoffset;
                    }
                }
            }
            // a[hint + lastoffset] <= key < a[hint + offset]
            offset = minimum(offset, maxoffset);

            lastoffset += hint;
            offset += hint;
        }
        a -= hint;
        LASSERT(-1<=lastoffset && lastoffset<offset);
        LASSERT(offset<=n);

        //a[lastoffset] <= key < a[offset]
        //binary search
        ++lastoffset;
        while(lastoffset<offset){
            s32 m = lastoffset + ((offset-lastoffset)>>1);
            if(compFunc_(key, a[m])){
                offset = m;
            } else{
                lastoffset = m+1;
            }
        }
        LASSERT(lastoffset == offset);
        return offset;
    }

    template<class T, class Comp>
    void Timsort<T, Comp>::merge_init(s32 size)
    {
        allocated_ = (size + 1)/2;
        if(MERGESTATE_TEMP_SIZE/2<allocated_){
            allocated_ = MERGESTATE_TEMP_SIZE/2;
        }
        a_.keys_ = temparray_;
        numPendings_ = 0;
        min_gallop_ = MIN_GALLOP;
    }

    template<class T, class Comp>
    void Timsort<T, Comp>::merge_free()
    {
        if(a_.keys_ != temparray_){
            LIME_FREE(a_.keys_);
        }
    }

    template<class T, class Comp>
    s32 Timsort<T, Comp>::merge_getmem(s32 need)
    {
        merge_free();

        a_.keys_ = (pointer_type)LIME_MALLOC(need * sizeof(value_type));
        if(NULL == a_.keys_){
            return -1;
        }
        allocated_ = need;
        return 0;
    }

    template<class T, class Comp>
    s32 Timsort<T, Comp>::merge_compute_minrun(s32 n)
    {
        LASSERT(0<=n);
        s32 r=0;
        while(64<=n){
            r |= n&1;
            n >>= 1;
        }
        return n+r;
    }

    template<class T, class Comp>
    void Timsort<T, Comp>::merge_collapse()
    {
        struct s_slice* p = pending_;
        while(1<numPendings_){
            s32 n = numPendings_ - 2;
            if(0<n && p[n-1].length_<=p[n].length_ + p[n+1].length_){
                if(p[n-1].length_<p[n+1].length_){
                    --n;
                }
                merge_at(n);

            } else if(p[n].length_<=p[n+1].length_){
                merge_at(n);

            } else {
                break;
            }
        }
    }

    template<class T, class Comp>
    void Timsort<T, Comp>::merge_force_collapse()
    {
        s_slice* p = pending_;
        while(1<numPendings_){
            s32 n = numPendings_ - 2;
            if(0<n && p[n-1].length_<p[n+1].length_){
                --n;
            }
            merge_at(n);
        }
    }

    // ssaとssbをマージ, na<=nb
    template<class T, class Comp>
    void Timsort<T, Comp>::merge_lo(sortslice ssa, s32 na, sortslice ssb, s32 nb)
    {
        LASSERT(0<na && 0<nb);
        LASSERT(ssa.keys_ + na == ssb.keys_);
        MERGE_GETMEM(na);

        //temporaryにコピー
        sortslice::memcpy(a_, 0, ssa, 0, na);
        sortslice dest = ssa;
        ssa = a_;

        //marge_atのgallop_rightで, ssb.keys_[0] < ssa.keys_[0]になっている
        sortslice::copy_incr(dest, ssb);
        --nb;

        if(0 == nb){
            goto Succeed;
        }
        if(1 == na){
            goto CopyB;
        }

        s32 min_gallop = min_gallop_;
        s32 k;
        for(;;){
            s32 acount = 0;
            s32 bcount = 0;

            //ssa, ssbどちらかから, 連続してmin_gallop個選択するまでマージ
            //値が同じ場合はssaを選択
            for(;;){
                LASSERT(1<na && 0<nb);
                if(compFunc_(ssb.keys_[0], ssa.keys_[0])){
                    sortslice::copy_incr(dest, ssb);
                    ++bcount;
                    acount = 0;
                    --nb;
                    if(0==nb){
                        goto Succeed;
                    }
                    if(min_gallop<=bcount){
                        break;
                    }
                } else{
                    sortslice::copy_incr(dest, ssa);
                    ++acount;
                    bcount = 0;
                    --na;
                    if(1 == na){
                        goto CopyB;
                    }
                    if(min_gallop<=acount){
                        break;
                    }
                } //if(compFunc_
            } //for(;;){

            ++min_gallop;
            do{
                LASSERT(1<na && 0<nb);
                min_gallop -= (1<min_gallop);
                min_gallop_ = min_gallop;

                //ssa.keys_[k-1]<= ssb.keys_[0] < ssa.keys_[k], 0<k
                k = gallop_right(ssb.keys_[0], ssa.keys_, na, 0);
                acount = k;
                if(k){
                    sortslice::memcpy(dest, 0, ssa, 0, k);
                    sortslice::advance(dest, k);
                    sortslice::advance(ssa, k);
                    na -= k;
                    if(1 == na){
                        goto CopyB;
                    }
                    //比較関数が一貫していればna==0にはならない
                    if(0 == na){
                        goto Succeed;
                    }
                }
                sortslice::copy_incr(dest, ssb);
                --nb;
                if(0 == nb){
                    goto Succeed;
                }
                //ssb.keys_[k-1]< ssa.keys_[0] <= ssb.keys_[k], 0<k
                k = gallop_left(ssa.keys_[0], ssb.keys_, nb, 0);
                bcount = k;
                if(k){
                    sortslice::memmove(dest, 0, ssb, 0, k);
                    sortslice::advance(dest, k);
                    sortslice::advance(ssb, k);
                    nb -= k;
                    if(0 == nb){
                        goto Succeed;
                    }
                }
                sortslice::copy_incr(dest, ssa);
                --na;
                if(1 == na){
                    goto CopyB;
                }
            } while(MIN_GALLOP <= acount || MIN_GALLOP <= bcount);
            ++min_gallop;
            min_gallop_ = min_gallop;
        }
Succeed:
        if(na){
            sortslice::memcpy(dest, 0, ssa, 0, na);
        }
        return;
CopyB:
        //ssbの残り, ssa[0]をコピー
        //ssbの残り < ssa[0]
        LASSERT(1 == na && 0<nb);
        sortslice::memmove(dest, 0, ssb, 0, nb);
        sortslice::copy(dest, nb, ssa, 0);
    }

    // ssaとssbをマージ, nb<=na
    template<class T, class Comp>
    void Timsort<T, Comp>::merge_hi(sortslice ssa, s32 na, sortslice ssb, s32 nb)
    {
        LASSERT(0<na && 0<nb);
        LASSERT(ssa.keys_ + na == ssb.keys_);

        MERGE_GETMEM(nb);

        //temporaryにコピー
        sortslice dest = ssb;
        sortslice::advance(dest, nb-1);
        sortslice::memcpy(a_, 0, ssb, 0, nb);
        sortslice basea = ssa;
        sortslice baseb = a_;
        ssb.keys_ = a_.keys_ + nb - 1;
        sortslice::advance(ssa, na-1);
        //marge_atのgallop_leftで, ssb.keys_[nb-1]< ssb.keys_[na-1] <= ssb.keys_[nb]
        sortslice::copy_decr(dest, ssa);
        --na;
        if(0 == na){
            goto Succeed;
        }
        if(1 == nb){
            goto CopyA;
        }

        s32 min_gallop = min_gallop_;
        s32 k;
        for(;;){
            s32 acount = 0;
            s32 bcount = 0;

            //ssa, ssbどちらかから, 連続してmin_gallop個選択するまでマージ
            //値が同じ場合はssaを選択
            for(;;){
                LASSERT(0<na && 1<nb);
                if(compFunc_(ssb.keys_[0], ssa.keys_[0])){
                    sortslice::copy_decr(dest, ssa);
                    ++acount;
                    bcount = 0;
                    --na;
                    if(0==na){
                        goto Succeed;
                    }
                    if(min_gallop<=acount){
                        break;
                    }
                } else{
                    sortslice::copy_decr(dest, ssb);
                    ++bcount;
                    acount = 0;
                    --nb;
                    if(1 == nb){
                        goto CopyA;
                    }
                    if(min_gallop<=bcount){
                        break;
                    }
                } //if(compFunc_
            } //for(;;){

            ++min_gallop;
            do{
                LASSERT(0<na && 1<nb);
                min_gallop -= (1<min_gallop);
                min_gallop_ = min_gallop;

                //basea.keys_[k-1]<= ssb.keys_[0] < basea.keys_[k], 0<k
                k = gallop_right(ssb.keys_[0], basea.keys_, na, na-1);
                k = na - k;
                acount = k;
                if(k){
                    sortslice::advance(dest, -k);
                    sortslice::advance(ssa, -k);
                    sortslice::memmove(dest, 1, ssa, 1, k);
                    na -= k;
                    if(0 == na){
                        goto Succeed;
                    }
                }
                sortslice::copy_decr(dest, ssb);
                --nb;
                if(1 == nb){
                    goto CopyA;
                }
                //baseb.keys_[k-1]< ssa.keys_[0] <= baseb.keys_[k], 0<k
                k = gallop_left(ssa.keys_[0], baseb.keys_, nb, nb-1);
                k = nb - k;
                bcount = k;
                if(k){
                    sortslice::advance(dest, -k);
                    sortslice::advance(ssb, -k);
                    sortslice::memcpy(dest, 1, ssb, 1, k);
                    nb -= k;
                    if(1 == nb){
                        goto CopyA;
                    }
                    //比較関数が一貫していればnb==0にはならない
                    if(0 == nb){
                        goto Succeed;
                    }
                }
                sortslice::copy_decr(dest, ssa);
                --na;
                if(0 == na){
                    goto Succeed;
                }
            } while(MIN_GALLOP <= acount || MIN_GALLOP <= bcount);
            ++min_gallop;
            min_gallop_ = min_gallop;
        }
Succeed:
        if(nb){
            sortslice::memcpy(dest, -(nb-1), baseb, 0, nb);
        }
        return;
CopyA:
        LASSERT(1 == nb && 0<na);
        sortslice::memmove(dest, 1-na, ssa, 1-na, na);
        sortslice::advance(dest, -na);
        sortslice::advance(ssa, -na);
        sortslice::copy(dest, 0, ssb, 0);
    }

    // スタックのn, n+1番目の部分列をマージ
    template<class T, class Comp>
    void Timsort<T, Comp>::merge_at(s32 n)
    {
        LASSERT(2<=numPendings_);
        LASSERT(0<=n);
        LASSERT(n == numPendings_-2 || n == numPendings_-3);
        LASSERT(0<pending_[n].length_ && 0<pending_[n+1].length_);

        sortslice ssa = pending_[n].base_;
        s32 na = pending_[n].length_;
        sortslice ssb = pending_[n+1].base_;
        s32 nb = pending_[n+1].length_;
        LASSERT(0<na && 0<nb);
        LASSERT(ssa.keys_ + na == ssb.keys_);

        //スタックn番目にマージ結果を保存
        //n == numPendings_-3 ならば, スタックのn+2をシフト
        if(n == numPendings_-3){
            pending_[n+1] = pending_[n+2];
        }
        pending_[n].length_ = na + nb;
        --numPendings_;

        s32 k;
        //ssa.keys_[k-1]<= ssb.keys_[0] < ssa.keys_[k]
        k = gallop_right(ssb.keys_[0], ssa.keys_, na, 0);
        sortslice::advance(ssa, k); //ssaのk-1番目までは, ssb以下
        na -= k;
        if(na == 0){
            return;
        }

        //ssb.keys_[nb-1]< ssb.keys_[na-1] <= ssb.keys_[nb]
        nb = gallop_left(ssa.keys_[na-1], ssb.keys_, nb, nb-1);
        if(nb<=0){
            return;
        }
        if(na<=nb){
            merge_lo(ssa, na, ssb, nb);
        } else{
            merge_hi(ssa, na, ssb, nb);
        }
    }

    template<class T, class Comp>
    void Timsort<T, Comp>::sort(pointer_type keys, s32 size)
    {
        LASSERT(NULL != keys);
        LASSERT(0<=size);

        if(size<2){
            return;
        }
        sortslice lo;
        lo.keys_ = keys;
        merge_init(size);
        s32 minrun = merge_compute_minrun(size);

        do{
            s32 descending;
            s32 n = count_run(lo.keys_, lo.keys_+size, descending);
            if(descending){
                sortslice::reverse(lo, n);
            }

            if(n<minrun){
                const s32 force = (size<=minrun)? size : minrun;
                binarysort(lo, lo.keys_ + force, lo.keys_+n);
                n = force;
            }
            LASSERT(numPendings_<MAX_MERGE_PENDING);
            pending_[numPendings_].base_ = lo;
            pending_[numPendings_].length_ = n;
            ++numPendings_;
            merge_collapse();
            sortslice::advance(lo, n);
            size -= n;
        } while(0<size);
        merge_force_collapse();
    }
}
#endif //INC_LRENDER_SORT_H__
