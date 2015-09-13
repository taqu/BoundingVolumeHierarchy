#ifndef INC_LRENDER_LRENDER_H__
#define INC_LRENDER_LRENDER_H__
/**
@file lrender.h
@author t-sakai
@date 2015/09/04 create
*/
#include <assert.h>
#include <malloc.h>
#include <new>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>

//#include <mmintrin.h>  //MMX命令セット
#include <xmmintrin.h> //SSE命令セット
#include <emmintrin.h> //SSE2命令セット

#define NOMINMAX
#include <Windows.h>
#include <sstream>

/// 16バイトアライメント変数指定
#ifdef _MSC_VER
#define LIME_ALIGN16 __declspec(align(16))
#define LIME_ALIGN(x) __declspec(align(x))
#else
#define LIME_ALIGN16 __attribute__((align(16)))
#define LIME_ALIGN(x) __attribute__((align(x)))
#endif

#define LIME_PLACEMENT_NEW(ptr) new(ptr)
#define LIME_DELETE(p) delete p; (p)=NULL
#define LIME_DELETE_NONULL(p) delete p

#define LIME_DELETE_ARRAY(p) delete[] (p); (p)=NULL

#define LIME_MALLOC(size) (::malloc(size))
#define LIME_FREE(mem) ::free(mem); (mem)=NULL

#define LIME_ALIGNED_MALLOC(size, align) (_aligned_malloc(size, align))
#define LIME_ALIGNED_FREE(mem, align) _aligned_free(mem); (mem)=NULL

// Assertion
//-------------------
#if defined(_DEBUG)

#if defined(ANDROID)
#define LASSERT(expression) {if((expression)==false){__android_log_assert("assert", "lime", "%s (%d)", __FILE__, __LINE__);}}while(0)

#elif defined(__GNUC__)
#define LASSERT(expression) ( assert(expression) )

#else
#define LASSERT(expression) ( assert(expression) )
#endif

#else
#define LASSERT(expression)
#endif

namespace lmath
{
    typedef __m128 lm128; /// XMMレジスタに対応した単精度浮動小数点型
    typedef __m128i lm128i;
    typedef __m64 lm64;
}

namespace lrender
{
    typedef char Char;
    typedef wchar_t WChar;
    typedef __int8 s8;
    typedef __int16 s16;
    typedef __int32 s32;
    typedef __int64 s64;

    typedef unsigned __int8 u8;
    typedef unsigned __int16 u16;
    typedef unsigned __int32 u32;
    typedef unsigned __int64 u64;

    typedef float f32;
    typedef double f64;

    typedef intptr_t  intptr_t;
    typedef uintptr_t  uintptr_t;
    typedef ptrdiff_t  ptrdiff_t;
    typedef size_t lsize_t;

    extern const f32 Epsilon;

    template<class T>
    inline T absolute(const T x)
    {
        return abs(x);
    }

    template<class T>
    void swap(T& x0, T& x1)
    {
        T tmp = x0;
        x0 = x1;
        x1 = tmp;
    }

    template<class T>
    T minimum(const T& x0, const T& x1)
    {
        return (x0<x1)? x0 : x1;
    }

    template<class T>
    T maximum(const T& x0, const T& x1)
    {
        return (x0<x1)? x1 : x0;
    }

    inline bool isZero(f32 x)
    {
        return (absolute(x) < Epsilon);
    }

    inline bool isZeroPositive(f32 x)
    {
        return (x < Epsilon);
    }

    inline bool isZeroNegative(f32 x)
    {
        return (-Epsilon < x);
    }

    template<class T>
    inline T clamp(T val, T low, T high)
    {
        if (val <= low) return low;
        else if (val >= high) return high;
        else return val;
    }

    inline f32 clamp01(f32 v)
    {
        s32* t = (s32*)&v;
        s32 s = (*t) >> 31;
        s = ~s;
        *t &= s;

        v -= 1.0f;
        s = (*t) >> 31;
        *t &= s;
        v += 1.0f;
        return v;
    }

    struct DefaultAllocator
    {
        static inline void* malloc(u32 size)
        {
            return LIME_MALLOC(size);
        }

        static inline void* malloc(u32 size, const char* /*file*/, int /*line*/)
        {
            return LIME_MALLOC(size);
        }

        static inline void* malloc(u32 size, u32 alignment, const char* /*file*/, int /*line*/)
        {
            return LIME_ALIGNED_MALLOC(size, alignment);
        }

        static inline void free(void* mem)
        {
            LIME_FREE(mem);
        }

        static inline void* malloc(u32 size, u32 alignment)
        {
            return LIME_ALIGNED_MALLOC(size, alignment);
        }

        static inline void free(void* mem, u32 alignment)
        {
            LIME_ALIGNED_FREE(mem, alignment);
        }
    };

    struct Align16Allocator
    {
        static inline void* malloc(u32 size)
        {
            return LIME_ALIGNED_MALLOC(size, 16);
        }

        static inline void* malloc(u32 size, const char* /*file*/, int /*line*/)
        {
            return LIME_ALIGNED_MALLOC(size, 16);
        }

        static inline void* malloc(u32 size, u32 alignment, const char* /*file*/, int /*line*/)
        {
            return LIME_ALIGNED_MALLOC(size, alignment);
        }

        static inline void free(void* mem)
        {
            LIME_ALIGNED_FREE(mem, 16);
        }

        static inline void* malloc(u32 size, u32 alignment)
        {
            return LIME_ALIGNED_MALLOC(size, alignment);
        }

        static inline void free(void* mem, u32 alignment)
        {
            LIME_ALIGNED_FREE(mem, alignment);
        }
    };


    extern const f32 HitEpsilon;

    struct RGB
    {
        f32 r_;
        f32 g_;
        f32 b_;
        f32 x_;
    };

    void printImage(const char* filename, RGB* rgb, s32 width, s32 height);

    class Vector3
    {
    public:
        Vector3()
        {}

        Vector3(f32 x, f32 y, f32 z)
            :x_(x)
            ,y_(y)
            ,z_(z)
        {}

        void zero()
        {
            x_ = y_ = z_ = 0.0f;
        }

        f32 operator[](s32 index) const
        {
            return reinterpret_cast<const f32*>(this)[index];
        }

        f32& operator[](s32 index)
        {
            return reinterpret_cast<f32*>(this)[index];
        }


        f32 length() const
        {
            return ::sqrtf(x_*x_ + y_*y_ + z_*z_);
        }

        f32 halfArea() const
        {
            return x_*y_ + y_*z_ + z_*x_;
        }

        friend static Vector3 operator+(const Vector3& v0, const Vector3& v1)
        {
            return Vector3(v0.x_+v1.x_, v0.y_+v1.y_, v0.z_+v1.z_);
        }

        friend static Vector3 operator-(const Vector3& v0, const Vector3& v1)
        {
            return Vector3(v0.x_-v1.x_, v0.y_-v1.y_, v0.z_-v1.z_);
        }

        friend static Vector3 operator*(const Vector3& v, f32 a)
        {
            return Vector3(a*v.x_, a*v.y_, a*v.z_);
        }

        friend static Vector3 operator*(f32 a, const Vector3& v)
        {
            return Vector3(a*v.x_, a*v.y_, a*v.z_);
        }

        f32 x_, y_, z_;
    };

    Vector3 normalize(const Vector3& v);
    f32 dot(const Vector3& v0, const Vector3& v1);
    Vector3 cross(const Vector3& v0, const Vector3& v1);

    class Vector4
    {
    public:

        f32 x_, y_, z_, w_;
    };

    class Ray;

    struct BBox
    {
        BBox()
        {}

        BBox(const Vector3& bmin, const Vector3& bmax)
            :bmin_(bmin)
            ,bmax_(bmax)
        {}

        void zero()
        {
            bmin_.zero();
            bmax_.zero();
        }

        void setInvalid()
        {
            bmin_.x_ = bmin_.y_ = bmin_.z_ = FLT_MAX;
            bmax_.x_ = bmax_.y_ = bmax_.z_ = -FLT_MAX;
        }

        Vector3 extent() const
        {
            return bmax_ - bmin_;
        }

        Vector3 diagonal() const
        {
            return Vector3(
                bmax_.x_-bmin_.x_,
                bmax_.y_-bmin_.y_,
                bmax_.z_-bmin_.z_);
        }

        void extend(const BBox& bbox);
        s32 maxExtentAxis() const;

        f32 halfArea() const;

        bool testRay(f32& tmin, f32& tmax, const Ray& ray) const;

        Vector3 bmin_;
        Vector3 bmax_;
    };

    enum PrimitiveType
    {
        Primitive_Point,
        Primitive_Face,
        Primitive_Sphere,
    };

    class Point
    {
    public:
        static const s32 Type = Primitive_Point;

        f32 getCentroidX() const;
        f32 getCentroidY() const;
        f32 getCentroidZ() const;
        Vector3 getCentroid() const;

        BBox getBBox() const;

        Vector3 position_;
    };

    class Face
    {
    public:
        static const s32 Type = Primitive_Face;

        f32 getCentroidX() const;
        f32 getCentroidY() const;
        f32 getCentroidZ() const;
        Vector3 getCentroid() const;

        BBox getBBox() const;

        Vector3 point_[3];
    };

    class Sphere
    {
    public:
        static const s32 Type = Primitive_Sphere;

        Sphere()
        {}

        Sphere(const Vector3& position, f32 radius)
            :position_(position)
            ,radius_(radius)
        {}

        f32 getCentroidX() const;
        f32 getCentroidY() const;
        f32 getCentroidZ() const;
        Vector3 getCentroid() const;

        BBox getBBox() const;

        bool testRay(f32& t, const Ray& ray) const;

        Vector3 position_;
        f32 radius_;
    };

    class Ray
    {
    public:
        Ray()
        {}

        Ray(const Vector3& origin,
            const Vector3& direction,
            f32 t);

        void invertDirection();

        void setDirection(const Vector3& direction);
        void setDirection(const Vector3& direction, const Vector3& invDirection);

        Vector3 origin_;
        Vector3 direction_;
        Vector3 invDirection_;
        f32 t_;
    };

    struct HitRecord
    {
        //void calcPoint(Vector3& point);
        //void calcNormal(Vector3& normal);

        f32 t_; //光線の半直線距離パラメータ
        f32 u_;
        f32 v_;
        s32 type_;
        const void* primitive_;
    };

    class vector_arena_dynamic_inc_size
    {
    public:
        static const s32 DEFAULT_INCREMENT_SIZE = 16;

        void setIncSize(s32 size)
        {
            LASSERT(0<size);
            incSize_ = size;
        }

    protected:
        vector_arena_dynamic_inc_size()
            :incSize_(DEFAULT_INCREMENT_SIZE)
        {}

        vector_arena_dynamic_inc_size(const vector_arena_dynamic_inc_size& rhs)
            :incSize_(rhs.incSize_)
        {}

        explicit vector_arena_dynamic_inc_size(s32 incSize)
            :incSize_(incSize)
        {}

        void swap(vector_arena_dynamic_inc_size& rhs)
        {
            lrender::swap(incSize_, rhs.incSize_);
        }

        s32 incSize_;
    };

    template<s32 INC_SIZE>
    class vector_arena_static_inc_size
    {
    public:
        static const s32 incSize_ = INC_SIZE;

        void setIncSize(s32 /*size*/)
        {
        }

    protected:
        vector_arena_static_inc_size()
        {}

        vector_arena_static_inc_size(const vector_arena_static_inc_size& /*rhs*/)
        {}

        explicit vector_arena_static_inc_size(s32 /*incSize*/)
        {}

        void swap(vector_arena_static_inc_size& /*rhs*/)
        {
        }
    };

    template<class T, class Allocator=DefaultAllocator, class IncSize=vector_arena_static_inc_size<16> >
    class vector_arena : public IncSize
    {
    public:
        typedef vector_arena<T, Allocator, IncSize> this_type;
        typedef s32 size_type;
        typedef T* iterator;
        typedef const T* const_iterator;
        typedef Allocator allocator_type;
        typedef IncSize inc_size_type;

        vector_arena();
        vector_arena(const this_type& rhs);
        explicit vector_arena(s32 incSize);
        vector_arena(s32 size, s32 incSize);
        ~vector_arena();

        s32 size() const{ return size_;}
        s32 capacity() const{ return capacity_;}

        T& operator[](s32 index)
        {
            LASSERT(0<=index && index<size_);
            return items_[index];
        }

        const T& operator[](s32 index) const
        {
            LASSERT(0<=index && index<size_);
            return items_[index];
        }

        T& front()
        {
            LASSERT(0<size_);
            return items_[0];
        }

        const T& front() const
        {
            LASSERT(0<size_);
            return items_[0];
        }

        T& back()
        {
            LASSERT(0<size_);
            return items_[size_-1];
        }

        const T& back() const
        {
            LASSERT(0<size_);
            return items_[size_-1];
        }


        void push_back(const T& t);
        void pop_back();

        iterator begin(){ return items_;}
        const_iterator begin() const{ return items_;}

        iterator end(){ return items_ + size_;}
        const_iterator end() const{ return items_ + size_;}

        void clear();
        void swap(this_type& rhs);
        void reserve(s32 capacity);
        void resize(s32 size);

        void removeAt(s32 index);
        s32 find(const T& ptr) const;
    private:
        this_type& operator=(const this_type&);

        s32 capacity_;
        s32 size_;
        T *items_;
    };

    template<class T, class Allocator, class IncSize>
    vector_arena<T, Allocator, IncSize>::vector_arena()
        :capacity_(0)
        ,size_(0)
        ,items_(NULL)
    {
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T, Allocator, IncSize>::vector_arena(const this_type& rhs)
        :inc_size_type(rhs)
        ,capacity_(rhs.capacity_)
        ,size_(rhs.size_)
    {
        items_ = reinterpret_cast<T*>(allocator_type::malloc(capacity_*sizeof(T)));
        for(s32 i=0; i<size_; ++i){
            LIME_PLACEMENT_NEW(&items_[i]) T(rhs.items_[i]);
        }
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T, Allocator, IncSize>::vector_arena(s32 incSize)
        :inc_size_type(incSize)
        ,capacity_(0)
        ,size_(0)
        ,items_(NULL)
    {
        LASSERT(incSize>0);
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T, Allocator, IncSize>::vector_arena(s32 size, s32 incSize)
        :inc_size_type(incSize)
        ,capacity_( (size>incSize)? size : incSize )
        ,size_(size)
    {
        LASSERT(incSize>0);

        items_ = reinterpret_cast<T*>(allocator_type::malloc(capacity_*sizeof(T)));

        for(s32 i=0; i<size_; ++i){
            LIME_PLACEMENT_NEW(&items_[i]) T();
        }
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T, Allocator, IncSize>::~vector_arena()
    {
        for(s32 i=0; i<size_; ++i){
            items_[i].~T();
        }
        allocator_type::free(items_);
        items_ = NULL;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T, Allocator, IncSize>::push_back(const T& t)
    {
        if(size_ >= capacity_){
            //新しいバッファ確保
            s32 newCapacity = capacity_ + inc_size_type::incSize_;
            T *newItems = reinterpret_cast<T*>( allocator_type::malloc(newCapacity*sizeof(T)) );

            //コピーコンストラクタでコピー。古い要素のデストラクト
            for(s32 i=0; i<size_; ++i){
                LIME_PLACEMENT_NEW(&newItems[i]) T(items_[i]);
                items_[i].~T();
            }
            LIME_PLACEMENT_NEW(&newItems[size_]) T(t);

            //古いバッファ破棄
            allocator_type::free(items_);

            items_ = newItems;
            ++size_;
            capacity_ = newCapacity;

        }else{
            LIME_PLACEMENT_NEW(&items_[size_]) T(t);
            ++size_;
        }
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T, Allocator, IncSize>::pop_back()
    {
        --size_;
        items_[size_].~T();
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T, Allocator, IncSize>::clear()
    {
        for(s32 i=0; i<size_; ++i){
            items_[i].~T();
        }
        size_ = 0;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T, Allocator, IncSize>::swap(this_type& rhs)
    {
        inc_size_type::swap(rhs);
        lcore::swap(capacity_, rhs.capacity_);
        lcore::swap(size_, rhs.size_);
        lcore::swap(items_, rhs.items_);
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T, Allocator, IncSize>::reserve(s32 capacity)
    {
        if(capacity<=capacity_){
            return;
        }

        //新しいバッファ確保
        T* newItems = reinterpret_cast<T*>( allocator_type::malloc(capacity*sizeof(T)) );

        //コピーコンストラクタでコピー。古い要素のデストラクト
        for(s32 i=0; i<size_; ++i){
            LIME_PLACEMENT_NEW(&newItems[i]) T(items_[i]);
            items_[i].~T();
        }

        //古いバッファ破棄
        allocator_type::free(items_);

        items_ = newItems;
        capacity_ = capacity;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T, Allocator, IncSize>::resize(s32 size)
    {
        if(size < size_){
            //デストラクト
            for(s32 i=size; i<size_; ++i){
                items_[i].~T();
            }

        }else{
            reserve(size);
            for(s32 i=size_; i<size; ++i){
                LIME_PLACEMENT_NEW(&items_[i]) T;
            }
        }
        size_ = size;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T, Allocator, IncSize>::removeAt(s32 index)
    {
        LASSERT(0<=index && index<size_);
        for(s32 i=index+1; i<size_; ++i){
            items_[i-1] = items_[i];
        }
        --size_;
        items_[size_].~T();
    }

    template<class T, class Allocator, class IncSize>
    s32 vector_arena<T, Allocator, IncSize>::find(const T& ptr) const
    {
        for(s32 i=0; i<size_; ++i){
            if(ptr == items_[i]){
                return i;
            }
        }
        return -1;
    }

    //-----------------------------------------------------------------
    //---
    //--- vector_arena ポインタ特殊化
    //---
    //-----------------------------------------------------------------
    template<class T, class Allocator, class IncSize>
    class vector_arena<T*, Allocator, IncSize> : public IncSize
    {
    public:
        typedef vector_arena<T, Allocator, IncSize> this_type;
        typedef s32 size_type;
        typedef T** iterator;
        typedef const T** const_iterator;
        typedef Allocator allocator_type;
        typedef IncSize inc_size_type;

        vector_arena();
        vector_arena(const this_type& rhs);
        explicit vector_arena(s32 incSize);
        vector_arena(s32 size, s32 incSize);
        ~vector_arena();

        s32 size() const{ return size_;}
        s32 capacity() const{ return capacity_;}

        T*& operator[](s32 index)
        {
            LASSERT(0<=index && index<size_);
            return items_[index];
        }

        const T* operator[](s32 index) const
        {
            LASSERT(0<=index && index<size_);
            return items_[index];
        }

        T*& front()
        {
            LASSERT(0<size_);
            return items_[0];
        }

        const T*& front() const
        {
            LASSERT(0<size_);
            return items_[0];
        }

        T*& back()
        {
            LASSERT(0<size_);
            return items_[size_-1];
        }

        const T*& back() const
        {
            LASSERT(0<size_);
            return items_[size_-1];
        }


        void push_back(T* const t);
        void pop_back();

        iterator begin(){ return items_;}
        const_iterator begin() const{ return items_;}

        iterator end(){ return items_ + size_;}
        const_iterator end() const{ return items_ + size_;}

        void clear();
        void swap(this_type& rhs);
        void reserve(s32 capacity);
        void resize(s32 size);

        void removeAt(s32 index);
        s32 find(const T* ptr) const;
    private:
        this_type& operator=(const this_type&);

        s32 capacity_;
        s32 size_;
        T** items_;
    };

    template<class T, class Allocator, class IncSize>
    vector_arena<T*, Allocator, IncSize>::vector_arena()
        :capacity_(0)
        ,size_(0)
        ,items_(NULL)
    {
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T*, Allocator, IncSize>::vector_arena(const this_type& rhs)
        :inc_size_type(rhs)
        ,capacity_(rhs.capacity_)
        ,size_(rhs.size_)
    {
        items_ = reinterpret_cast<T**>(allocator_type::malloc(capacity_*sizeof(T*)));

        for(s32 i=0; i<size_; ++i){
            items_[i] = NULL;
        }
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T*, Allocator, IncSize>::vector_arena(s32 incSize)
        :inc_size_type(incSize)
        ,capacity_(0)
        ,size_(0)
        ,items_(NULL)
    {
        LASSERT(incSize>0);
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T*, Allocator, IncSize>::vector_arena(s32 size, s32 incSize)
        :inc_size_type(incSize)
        ,capacity_( (size>incSize)? size : incSize )
        ,size_(size)
    {
        LASSERT(incSize>0);
        items_ = reinterpret_cast<T**>(allocator_type::malloc(capacity_*sizeof(T*)));

        for(s32 i=0; i<size_; ++i){
            items_[i] = NULL;
        }
    }

    template<class T, class Allocator, class IncSize>
    vector_arena<T*, Allocator, IncSize>::~vector_arena()
    {
        allocator_type::free(items_);
        items_ = NULL;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T*, Allocator, IncSize>::push_back(T* const t)
    {
        if(size_ >= capacity_){
            //新しいバッファ確保
            s32 newCapacity = capacity_ + inc_size_type::incSize_;
            T** newItems = reinterpret_cast<T**>(allocator_type::malloc(newCapacity*sizeof(T*)));

            //コピー
            for(s32 i=0; i<size_; ++i){
                newItems[i] = items_[i];
            }
            //古いバッファ破棄
            allocator_type::free(items_);

            items_ = newItems;
            capacity_ = newCapacity;

        }
        items_[size_] = t;
        ++size_;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T*, Allocator, IncSize>::pop_back()
    {
        LASSERT(size_>0);
        --size_;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T*, Allocator, IncSize>::clear()
    {
        size_ = 0;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T*, Allocator, IncSize>::swap(this_type& rhs)
    {
        inc_size_type::swap(rhs);
        lcore::swap(capacity_, rhs.capacity_);
        lcore::swap(size_, rhs.size_);
        lcore::swap(items_, rhs.items_);
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T*, Allocator, IncSize>::reserve(s32 capacity)
    {
        if(capacity<=capacity_){
            return;
        }

        //新しいバッファ確保
        T** newItems = reinterpret_cast<T**>(allocator_type::malloc(capacity*sizeof(T*)));

        //コピー
        for(s32 i=0; i<size_; ++i){
            newItems[i] = items_[i];
        }

        //古いバッファ破棄
        allocator_type::free(items_);

        items_ = newItems;
        capacity_ = capacity;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T*, Allocator, IncSize>::resize(s32 size)
    {
        if(size > size_){
            reserve(size);
        }
        size_ = size;
    }

    template<class T, class Allocator, class IncSize>
    void vector_arena<T*, Allocator, IncSize>::removeAt(s32 index)
    {
        LASSERT(0<=index && index<size_);
        for(s32 i=index+1; i<size_; ++i){
            items_[i-1] = items_[i];
        }
        --size_;
        items_[size_] = NULL;
    }

    template<class T, class Allocator, class IncSize>
    s32 vector_arena<T*, Allocator, IncSize>::find(const T* ptr) const
    {
        for(s32 i=0; i<size_; ++i){
            if(ptr == items_[i]){
                return i;
            }
        }
        return -1;
    }

#if defined(ANDROID) || defined(__GNUC__)
    typedef clock_t ClockType;
#else
    typedef u64 ClockType;
#endif

    void sleep(u32 milliSeconds);

    /// カウント取得
    ClockType getPerformanceCounter();

    /// 秒間カウント数
    ClockType getPerformanceFrequency();

    /// 秒単位の時間差分計算
    f64 calcTime64(ClockType prevTime, ClockType currentTime);

    inline f32 calcTime(ClockType prevTime, ClockType currentTime)
    {
        return static_cast<f32>(calcTime64(prevTime, currentTime));
    }

    /// ミリ秒単位の時間を取得
    u32 getTime();


    template<bool enable>
    struct Timer
    {
        Timer()
            :time_(0)
            ,count_(0)
            ,totalTime_(0.0f)
        {}

        void begin()
        {
            time_ = getPerformanceCounter();
        }

        void end()
        {
            totalTime_ += calcTime64(time_, getPerformanceCounter());
            ++count_;
        }

        f64 getAverage() const
        {
            return (0 == count_)? 0.0 : totalTime_/count_;
        }

        void reset();

        ClockType time_;
        s32 count_;
        f64 totalTime_;
    };

    template<bool enable>
    void Timer<enable>::reset()
    {
        time_ = 0;
        count_ = 0;
        totalTime_ = 0.0f;
    }

    template<>
    struct Timer<false>
    {
        void begin(){}
        void end(){}
        f64 getAverage() const{return 0.0;}
        void reset(){}
    };

    u32 getDefaultSeed();

    //---------------------------------------------
    //---
    //--- RandomXorshift
    //---
    //---------------------------------------------
    class RandomXorshift
    {
    public:
        RandomXorshift();
        explicit RandomXorshift(u32 seed);
        ~RandomXorshift();

        /**
        @brief 擬似乱数生成器初期化
        @param seed
        */
        void srand(u32 seed);

        /**
        @brief 0 - 0xFFFFFFFFUの乱数生成
        */
        u32 rand();

        /**
        @brief 0 - 1の乱数生成
        */
        f32 frand();

        /**
        @brief 0.0 - 0.999999881の乱数生成
        */
        f32 frand2();

        /**
        @brief 0 - 1の乱数生成
        */
        f64 drand();

    private:
        u32 rand(u32 v, u32 i);

        u32 x_;
        u32 y_;
        u32 z_;
        u32 w_;
    };
}


#include "Sort.h"

namespace lrender
{
    struct SortFuncCentroid
    {
        SortFuncCentroid(const f32* centroids)
            :centroids_(centroids)
        {}

        bool operator()(s32 i0, s32 i1) const
        {
            return centroids_[i0] < centroids_[i1];
        }

        const f32* centroids_;
    };

    template<class T>
    class PrimitivePolicy
    {
    public:

        static inline Vector3 getCentroid(const T& primitive)
        {
            return primitive.getCentroid();
        }

        static inline BBox getBBox(const T& primitive)
        {
            return primitive.getBBox();
        }

        static inline bool testRay(f32& t, const T& primitive, const Ray& ray)
        {
            return primitive.testRay(t, ray);
        }

        static inline void sort(s32 numPrimitives, s32* primitiveIndices, const f32* centroids)
        {
            lrender::introsort(numPrimitives, primitiveIndices, SortFuncCentroid(centroids));
        }

        static inline void insertionsort(s32 numPrimitives, s32* primitiveIndices, const f32* centroids)
        {
            lrender::insertionsort(numPrimitives, primitiveIndices, SortFuncCentroid(centroids));
        }
    };

namespace qbvh
{
    //-----------------------------------------------------------
    // 線分とAABBの交差判定
    s32 testRayAABB(
        __m128 tmin,
        __m128 tmax,
        __m128 origin[3],
        __m128 invDir[3],
        const s32 sign[3],
        const __m128 bbox[2][3]);

    //AABBの交差判定
    s32 testAABB(const __m128 bbox0[2][3], const __m128 bbox1[2][3]);

    //
    s32 testSphereAABB(const __m128 position[3], const __m128 radius, const __m128 bbox[2][3]);
}
}
#endif //INC_LRENDER_LRENDER_H__
