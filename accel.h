#ifndef INC_ACCEL_ACCEL_H__
#define INC_ACCEL_ACCEL_H__
/**
@file accel.h
@author t-sakai
@date 2018/01/22 create
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

#include <utility>

//C++98 199711L
//C++11 201103L
#ifdef __cplusplus
#   if 201103L<=__cplusplus || 1900<=_MSC_VER
#       define ACC_CPP11 1
#   endif
#endif

#ifdef NULL
#   ifdef __cplusplus
#       ifdef ACC_CPP11
#           undef NULL
#           define NULL nullptr
#       endif
#   endif
#else
#   ifdef __cplusplus
#       ifdef ACC_CPP11
#           define NULL nullptr
#       else
#           define NULL 0
#       endif
#   else
#       define NULL ((void*)0)
#   endif
#endif

/// 16バイトアライメント変数指定
#ifdef _MSC_VER
#define ACC_ALIGN16 __declspec(align(16))
#define ACC_ALIGN(x) __declspec(align(x))
#else
#define ACC_ALIGN16 __attribute__((align(16)))
#define ACC_ALIGN(x) __attribute__((align(x)))
#endif

#define ACC_PLACEMENT_NEW(ptr) new(ptr)
#define ACC_DELETE(p) delete p; (p)=NULL
#define ACC_DELETE_NONULL(p) delete p

#define ACC_DELETE_ARRAY(p) delete[] (p); (p)=NULL

#define ACC_MALLOC(size) (::malloc(size))
#define ACC_FREE(mem) ::free(mem); (mem)=NULL

#define ACC_ALIGNED_MALLOC(size, align) (_aligned_malloc(size, align))
#define ACC_ALIGNED_FREE(mem, align) _aligned_free(mem); (mem)=NULL

// Assertion
//-------------------
#if defined(_DEBUG)

#if defined(ANDROID)
#define ACC_ASSERT(expression) {if((expression)==false){__android_log_assert("assert", "lime", "%s (%d)", __FILE__, __LINE__);}}while(0)

#elif defined(__GNUC__)
#define ACC_ASSERT(expression) ( assert(expression) )

#else
#define ACC_ASSERT(expression) ( assert(expression) )
#endif

#else
#define ACC_ASSERT(expression)
#endif

namespace accel
{
    typedef char Char;
    typedef unsigned char UChar;
    typedef char16_t Char16;
    typedef wchar_t WChar;
    typedef int8_t s8;
    typedef int16_t s16;
    typedef int32_t s32;
    typedef int64_t s64;

    typedef uint8_t u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    typedef float f32;
    typedef double f64;

    typedef intptr_t  intptr_t;
    typedef uintptr_t  uintptr_t;
    typedef ptrdiff_t  ptrdiff_t;
    typedef size_t size_t;

    typedef __m128 lm128; /// XMMレジスタに対応した単精度浮動小数点型
    typedef __m128i lm128i;
    typedef __m64 lm64;

#if defined(ANDROID)
    static constexpr f32 F32_EPSILON = 1.192092896e-07F;
    static constexpr f32 F64_EPSILON = 2.2204460492503131e-016;
#else
    static constexpr f32 F32_EPSILON = FLT_EPSILON;
    static constexpr f32 F64_EPSILON = DBL_EPSILON;
#endif

    static constexpr f32 F32_HITEPSILON = 1.0e-5f;

    enum Result
    {
        Result_Fail = 0,
        Result_Success = (0x01U<<0),
        Result_Front = (0x01U<<0)|(0x01U<<1),
        Result_Back = (0x01U<<0)|(0x01U<<2),
    }; 

    using std::move;

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

    inline bool isEqual(f32 x0, f32 x1)
    {
        return (absolute(x0-x1) < F32_EPSILON);
    }

    inline bool isZero(f32 x)
    {
        return (absolute(x) < F32_EPSILON);
    }

    inline bool isZeroPositive(f32 x)
    {
        return (x < F32_EPSILON);
    }

    inline bool isZeroNegative(f32 x)
    {
        return (-F32_EPSILON < x);
    }

    template<class T>
    inline T clamp(T val, T low, T high)
    {
        if (val <= low) return low;
        else if (val >= high) return high;
        else return val;
    }

    f32 clamp01(f32 v);

    // Returned value is undefined, if x==0
    u32 leadingzero(u32 x);

    struct RGB
    {
        f32 r_;
        f32 g_;
        f32 b_;
        f32 x_;
    };

    void printImage(const char* filename, RGB* rgb, s32 width, s32 height);

    //--- Vector3
    //--------------------------------------------------------------
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

        Vector3& operator*=(f32 a)
        {
            x_ *= a;
            y_ *= a;
            z_ *= a;
            return *this;
        }

        Vector3& operator*=(const Vector3& v)
        {
            x_ *= v.x_;
            y_ *= v.y_;
            z_ *= v.z_;
            return *this;
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

    //--- Vector4
    //--------------------------------------------------------------
    class Vector4
    {
    public:

        f32 x_, y_, z_, w_;
    };

    //--- Morton Code
    //--------------------------------------------------------------
    /**
    @brief Generate morton code 10 bits for each axis
    */
    s32 mortonCode3(s32 x, s32 y, s32 z);
    /**
    @brief Reconstruct each axis' values from morton code
    */
    void rmortonCode3(s32& x, s32& y, s32& z, s32 w);

    class Ray;

    //--- AABB
    //--------------------------------------------------------------
    struct AABB
    {
        AABB()
        {}

        AABB(const Vector3& bmin, const Vector3& bmax)
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

        void extend(const AABB& bbox);
        s32 maxExtentAxis() const;

        f32 halfArea() const;

        bool testRay(f32& tmin, f32& tmax, const Ray& ray) const;

        Vector3 bmin_;
        Vector3 bmax_;
    };

#ifdef _DEBUG
    inline bool operator==(const AABB& b0, const AABB& b1)
    {
        for(int i=0; i<3; ++i){
            if(!isEqual(b0.bmin_[i], b1.bmin_[i])){
                return false;
            }
            if(!isEqual(b0.bmax_[i], b1.bmax_[i])){
                return false;
            }
        }
        return true;
    }
#endif

    //--- PrimitiveType
    //--------------------------------------------------------------
    enum PrimitiveType
    {
        Primitive_Point,
        Primitive_Face,
        Primitive_Sphere,
    };

    //--- Point
    //--------------------------------------------------------------
    class Point
    {
    public:
        static const s32 Type = Primitive_Point;

        f32 getCentroidX() const;
        f32 getCentroidY() const;
        f32 getCentroidZ() const;
        Vector3 getCentroid() const;
        f32 getCentroid(s32 axis) const;

        AABB getBBox() const;

        Vector3 position_;
    };

    //--- Face
    //--------------------------------------------------------------
    class Face
    {
    public:
        static const s32 Type = Primitive_Face;

        f32 getCentroidX() const;
        f32 getCentroidY() const;
        f32 getCentroidZ() const;
        Vector3 getCentroid() const;
        f32 getCentroid(s32 axis) const;

        AABB getBBox() const;

        bool testRay(f32& t, const Ray& ray) const;
        Vector3 getNormal() const;

        Vector3 point_[3];
    };

    //--- Sphere
    //--------------------------------------------------------------
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
        f32 getCentroid(s32 axis) const;

        AABB getBBox() const;

        bool testRay(f32& t, const Ray& ray) const;

        Vector3 position_;
        f32 radius_;
    };

    //--- Ray
    //--------------------------------------------------------------
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

    //--- HitRecord
    //--------------------------------------------------------------
    struct HitRecord
    {
        f32 t_; //光線の半直線距離パラメータ
        const void* primitive_;
    };

    template<s32 size=16>
    struct ArrayStaticCapacityIncrement
    {
        static const s32 IncrementSize = size;

        static s32 getInitCapacity(s32 capacity)
        {
            return (capacity + (IncrementSize-1)) & ~(IncrementSize-1);
        }

        static s32 getNewCapacity(s32 capacity)
        {
            return capacity + IncrementSize;
        }
    };

    //-------------------------------------------------------
    //---
    //---
    //---
    //-------------------------------------------------------
    template<class T, class CapacityIncrement, class isTriviallyCopyable>
    class ArrayBuffer
    {
    protected:
        typedef ArrayBuffer<T, CapacityIncrement, isTriviallyCopyable> this_type;
        typedef T value_type;
        typedef s32 size_type;
        typedef CapacityIncrement capacity_increment_type;

        ArrayBuffer(const this_type&) = delete;
        this_type& operator=(const this_type&) = delete;

        ArrayBuffer()
            :capacity_(0)
            ,size_(0)
            ,items_(NULL)
        {}

        ArrayBuffer(this_type&& rhs)
            :capacity_(rhs.capacity_)
            ,size_(rhs.size_)
            ,items_(rhs.items_)
        {
            rhs.capacity_ = 0;
            rhs.size_ = 0;
            rhs.items_ = NULL;
        }

        explicit ArrayBuffer(size_type capacity)
            :capacity_(capacity_increment_type::getInitCapacity(capacity))
            ,size_(0)
            ,items_(NULL)
        {
            ACC_ASSERT(0<=capacity_);
            items_ = reinterpret_cast<T*>(ACC_ALIGNED_MALLOC(capacity_*sizeof(T), 16));
        }

        ~ArrayBuffer()
        {
            for(s32 i = 0; i<size_; ++i){
                items_[i].~value_type();
            }
            ACC_ALIGNED_FREE(items_, 16);
        }

        void helper_push_back(const value_type& t)
        {
            if(capacity_<=size_){
                helper_reserve(capacity_increment_type::getNewCapacity(capacity_));
            }
            ACC_PLACEMENT_NEW(&items_[size_]) value_type(t);
            ++size_;
        }
        void helper_push_back(value_type&& t)
        {
            if(capacity_<=size_){
                helper_reserve(capacity_increment_type::getNewCapacity(capacity_));
            }
            ACC_PLACEMENT_NEW(&items_[size_]) value_type(move(t));
            ++size_;
        }

        void helper_pop_back()
        {
            ACC_ASSERT(0<size_);
            --size_;
            items_[size_].~value_type();
        }

        void helper_clear()
        {
            for(s32 i = 0; i<size_; ++i){
                items_[i].~value_type();
            }
            size_ = 0;
        }

        void helper_reserve(size_type capacity)
        {
            if(capacity<=capacity_){
                return;
            }

            capacity = capacity_increment_type::getInitCapacity(capacity);
            value_type* newItems = reinterpret_cast<value_type*>(ACC_ALIGNED_MALLOC(capacity*sizeof(value_type), 16));

            //Copy, and destruct
            for(s32 i = 0; i<size_; ++i){
                ACC_PLACEMENT_NEW(&newItems[i]) value_type(move(items_[i]));
                items_[i].~value_type();
            }

            ACC_ALIGNED_FREE(items_, 16);

            items_ = newItems;
            capacity_ = capacity;
        }

        void helper_resize(size_type size)
        {
            if(size < size_){
                //Destruct redundant items
                for(s32 i = size; i<size_; ++i){
                    items_[i].~value_type();
                }

            } else{
                //Construct new items by default constructor
                helper_reserve(size);
                for(s32 i = size_; i<size; ++i){
                    ACC_PLACEMENT_NEW(&items_[i]) value_type;
                }
            }
            size_ = size;
        }

        void helper_removeAt(s32 index)
        {
            ACC_ASSERT(0<=index && index<size_);
            for(s32 i = index+1; i<size_; ++i){
                items_[i-1] = move(items_[i]);
            }
            --size_;
            items_[size_].~value_type();
        }

        size_type capacity_;
        size_type size_;
        value_type* items_;
    };


    template<class T, class CapacityIncrement>
    class ArrayBuffer<T, CapacityIncrement, std::true_type>
    {
    protected:
        typedef ArrayBuffer<T, CapacityIncrement, std::true_type> this_type;
        typedef T value_type;
        typedef s32 size_type;
        typedef CapacityIncrement capacity_increment_type;

        ArrayBuffer(const this_type&) = delete;
        this_type& operator=(const this_type&) = delete;

        ArrayBuffer()
            :capacity_(0)
            ,size_(0)
            ,items_(NULL)
        {}

        ArrayBuffer(this_type&& rhs)
            :capacity_(rhs.capacity_)
            ,size_(rhs.size_)
            ,items_(rhs.items_)
        {
            rhs.capacity_ = 0;
            rhs.size_ = 0;
            rhs.items_ = NULL;
        }

        explicit ArrayBuffer(size_type capacity)
            :capacity_(capacity_increment_type::getInitCapacity(capacity))
            ,size_(0)
            ,items_(NULL)
        {
            ACC_ASSERT(0<=capacity_);
            items_ = reinterpret_cast<T*>(ACC_ALIGNED_MALLOC(capacity_*sizeof(T), 16));
        }

        ~ArrayBuffer()
        {
            for(s32 i = 0; i<size_; ++i){
                items_[i].~value_type();
            }
            ACC_ALIGNED_FREE(items_, 16);
        }

        void helper_push_back(const value_type& t)
        {
            if(capacity_<=size_){
                helper_reserve(capacity_increment_type::getNewCapacity(capacity_));
            }
            ACC_PLACEMENT_NEW(&items_[size_]) value_type(t);
            ++size_;
        }
        void helper_push_back(value_type&& t)
        {
            if(capacity_<=size_){
                helper_reserve(capacity_increment_type::getNewCapacity(capacity_));
            }
            ACC_PLACEMENT_NEW(&items_[size_]) value_type(move(t));
            ++size_;
        }

        void helper_pop_back()
        {
            ACC_ASSERT(0<size_);
            --size_;
            items_[size_].~value_type();
        }

        void helper_clear()
        {
            for(s32 i = 0; i<size_; ++i){
                items_[i].~value_type();
            }
            size_ = 0;
        }

        void helper_reserve(size_type capacity)
        {
            if(capacity<=capacity_){
                return;
            }

            capacity = capacity_increment_type::getInitCapacity(capacity);
            value_type* newItems = reinterpret_cast<value_type*>(ACC_ALIGNED_MALLOC(capacity*sizeof(value_type), 16));

            //Copy, and destruct
            for(s32 i = 0; i<size_; ++i){
                ACC_PLACEMENT_NEW(&newItems[i]) value_type(items_[i]);
                items_[i].~value_type();
            }

            ACC_ALIGNED_FREE(items_, 16);

            items_ = newItems;
            capacity_ = capacity;
        }

        void helper_resize(size_type size)
        {
            if(size < size_){
                //Destruct redundant items
                for(s32 i = size; i<size_; ++i){
                    items_[i].~value_type();
                }

            } else{
                //Construct new items by default constructor
                helper_reserve(size);
                for(s32 i = size_; i<size; ++i){
                    ACC_PLACEMENT_NEW(&items_[i]) value_type;
                }
            }
            size_ = size;
        }

        void helper_removeAt(s32 index)
        {
            ACC_ASSERT(0<=index && index<size_);
            for(s32 i = index+1; i<size_; ++i){
                items_[i-1] = move(items_[i]);
            }
            --size_;
            items_[size_].~value_type();
        }

        size_type capacity_;
        size_type size_;
        value_type* items_;
    };

    //-------------------------------------------------------
    //---
    //---
    //---
    //-------------------------------------------------------
    template<class T, class CapacityIncrement=ArrayStaticCapacityIncrement<> >
    class Array : public ArrayBuffer<T, CapacityIncrement, typename std::is_trivially_copyable<T>::type>
    {
    public:
        typedef Array<T, CapacityIncrement> this_type;
        typedef ArrayBuffer<T, CapacityIncrement, typename std::is_trivially_copyable<T>::type> parent_type;
        typedef T value_type;
        typedef T* iterator;
        typedef const T* const_iterator;
        typedef s32 size_type;
        typedef CapacityIncrement capacity_increment_type;

        /**
        @return if lhs<rhs then true else false
        */
        typedef bool(*SortCmp)(const T& lhs, const T& rhs);

        Array();
        Array(this_type&& rhs);
        explicit Array(size_type capacity);
        ~Array();

        inline size_type capacity() const;
        inline size_type size() const;

        inline T& operator[](s32 index);
        inline const T& operator[](s32 index) const;

        inline T& front();
        inline const T& front() const;
        inline T& back();
        inline const T& back() const;

        void push_back(const T& t);
        void push_back(T&& t);
        void pop_back();

        inline iterator begin();
        inline const_iterator begin() const;

        inline iterator end();
        inline const_iterator end() const;

        void clear();
        void reserve(size_type capacity);
        void resize(size_type size);
        void removeAt(s32 index);
        void swap(this_type& rhs);

        this_type& operator=(this_type&& rhs);

        s32 find(const T& ptr) const;
        void insertionsort(const T& t, SortCmp cmp);
    private:
        Array(const this_type&) = delete;
        this_type& operator=(const this_type&) = delete;
    };

    template<class T, class CapacityIncrement>
    Array<T, CapacityIncrement>::Array()
    {
    }

    template<class T, class CapacityIncrement>
    Array<T, CapacityIncrement>::Array(this_type&& rhs)
        :parent_type(rhs)
    {
    }

    template<class T, class CapacityIncrement>
    Array<T, CapacityIncrement>::Array(size_type capacity)
        :parent_type(capacity)
    {
    }

    template<class T, class CapacityIncrement>
    Array<T, CapacityIncrement>::~Array()
    {
    }

    template<class T, class CapacityIncrement>
    inline typename Array<T, CapacityIncrement>::size_type
        Array<T, CapacityIncrement>::capacity() const
    {
        return capacity_;
    }

    template<class T, class CapacityIncrement>
    inline typename Array<T, CapacityIncrement>::size_type
        Array<T, CapacityIncrement>::size() const
    {
        return size_;
    }

    template<class T, class CapacityIncrement>
    inline T& Array<T, CapacityIncrement>::operator[](s32 index)
    {
        ACC_ASSERT(0<=index && index<size_);
        return items_[index];
    }

    template<class T, class CapacityIncrement>
    inline const T& Array<T, CapacityIncrement>::operator[](s32 index) const
    {
        ACC_ASSERT(0<=index && index<size_);
        return items_[index];
    }

    template<class T, class CapacityIncrement>
    inline T& Array<T, CapacityIncrement>::front()
    {
        ACC_ASSERT(0<size_);
        return items_[0];
    }

    template<class T, class CapacityIncrement>
    inline const T& Array<T, CapacityIncrement>::front() const
    {
        ACC_ASSERT(0<size_);
        return items_[0];
    }

    template<class T, class CapacityIncrement>
    inline T& Array<T, CapacityIncrement>::back()
    {
        ACC_ASSERT(0<size_);
        return items_[size_-1];
    }

    template<class T, class CapacityIncrement>
    inline const T& Array<T, CapacityIncrement>::back() const
    {
        ACC_ASSERT(0<size_);
        return items_[size_-1];
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::push_back(const T& t)
    {
        helper_push_back(t);
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::push_back(T&& t)
    {
        helper_push_back(move(t));
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::pop_back()
    {
        helper_pop_back();
    }

    template<class T, class CapacityIncrement>
    inline typename Array<T, CapacityIncrement>::iterator Array<T, CapacityIncrement>::begin()
    {
        return items_;
    }
    template<class T, class CapacityIncrement>
    inline typename Array<T, CapacityIncrement>::const_iterator Array<T, CapacityIncrement>::begin() const
    {
        return items_;
    }

    template<class T, class CapacityIncrement>
    inline typename Array<T, CapacityIncrement>::iterator Array<T, CapacityIncrement>::end()
    {
        return items_ + size_;
    }
    template<class T, class CapacityIncrement>
    inline typename Array<T, CapacityIncrement>::const_iterator Array<T, CapacityIncrement>::end() const
    {
        return items_ + size_;
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::clear()
    {
        helper_clear();
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::reserve(size_type capacity)
    {
        helper_reserve(capacity);
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::resize(size_type size)
    {
        helper_resize(size);
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::removeAt(s32 index)
    {
        helper_removeAt(index);
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::swap(this_type& rhs)
    {
        accel::swap(capacity_, rhs.capacity_);
        accel::swap(size_, rhs.size_);
        accel::swap(items_, rhs.items_);
    }

    template<class T, class CapacityIncrement>
    typename Array<T, CapacityIncrement>::this_type& Array<T, CapacityIncrement>::operator=(this_type&& rhs)
    {
        if(this == &rhs){
            return;
        }
        for(s32 i = 0; i<size_; ++i){
            items_[i].~value_type();
        }
        ACC_ALIGNED_FREE(items_, 16);

        capacity_ = rhs.capacity_;
        size_ = rhs.size_;
        items_ = rhs.items_;

        rhs.capacity_ = 0;
        rhs.size_ = 0;
        rhs.items_ = NULL;
    }

    template<class T, class CapacityIncrement>
    s32 Array<T, CapacityIncrement>::find(const T& ptr) const
    {
        for(s32 i=0; i<size_; ++i){
            if(ptr == items_[i]){
                return i;
            }
        }
        return -1;
    }

    template<class T, class CapacityIncrement>
    void Array<T, CapacityIncrement>::insertionsort(const T& t, SortCmp cmp)
    {
        s32 size = size_;
        push_back(t);
        for(s32 i=size-1; 0<=i; --i){
            if(cmp(items_[i+1], items_[i])){
                lray::swap(items_[i+1], items[i]);
                continue;
            }
            break;
        }
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

    template<bool enable>
    struct Timer
    {
        Timer()
            :time_(0)
            ,count_(0)
            ,totalTime_(0.0f)
        {}

        void restart()
        {
            reset();
            start();
        }

        void start()
        {
            time_ = getPerformanceCounter();
        }

        void stop()
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
        void restart(){}
        void start(){}
        void stop(){}
        f64 getAverage() const{return 0.0;}
        void reset(){}
    };

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

namespace accel
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

    struct SortFuncCode
    {
        SortFuncCode(const u32* codes)
            :codes_(codes)
        {}

        bool operator()(s32 i0, s32 i1) const
        {
            return codes_[i0] < codes_[i1];
        }

        const u32* codes_;
    };

    template<class T>
    class PrimitivePolicy
    {
    public:

        static inline Vector3 getCentroid(const T& primitive)
        {
            return primitive.getCentroid();
        }

        static inline f32 getCentroid(const T& primitive, s32 axis)
        {
            return primitive.getCentroid(axis);
        }

        static inline AABB getBBox(const T& primitive)
        {
            return primitive.getBBox();
        }

        static inline bool testRay(f32& t, const T& primitive, const Ray& ray)
        {
            return primitive.testRay(t, ray);
        }

        static inline void sort(s32 numPrimitives, s32* primitiveIndices, const f32* centroids)
        {
            introsort(numPrimitives, primitiveIndices, SortFuncCentroid(centroids));
        }

        static inline void insertionsort(s32 numPrimitives, s32* primitiveIndices, const f32* centroids)
        {
            accel::insertionsort(numPrimitives, primitiveIndices, SortFuncCentroid(centroids));
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
#endif //INC_ACCEL_ACCEL_H__
