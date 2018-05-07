/**
@file accel.cpp
@author t-sakai
@date 2018/01/22 create
*/
#include "accel.h"
#include <stdio.h>

namespace accel
{
    f32 clamp01(f32 v)
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

    u32 leadingzero(u32 x)
    {
#if defined(_MSC_VER)
        unsigned long n;
        _BitScanReverse(&n, x);
        return 31-n;
#elif defined(__GNUC__)
        return __builtin_clz(x);
#else
        u32 n = 0;
        if(x<=0x0000FFFFU){ n+=16; x<<=16;}
        if(x<=0x00FFFFFFU){ n+= 8; x<<= 8;}
        if(x<=0x0FFFFFFFU){ n+= 4; x<<= 4;}
        if(x<=0x3FFFFFFFU){ n+= 2; x<<= 2;}
        if(x<=0x7FFFFFFFU){ ++n;}
        return n;
#endif
    }


    void printImage(const char* filename, RGB* rgb, s32 width, s32 height)
    {
        ACC_ASSERT(NULL != filename);
        ACC_ASSERT(NULL != rgb);

        FILE* file = NULL;
        fopen_s(&file, filename, "wb");
        if(NULL == file){
            return;
        }

        fprintf(file, "P3\n");
        fprintf(file, "%d %d\n", width, height);
        fprintf(file, "255\n");

        for(s32 i=0; i<height; ++i){
            for(s32 j=0; j<width; ++j){
                const RGB& pixel = rgb[width*i + j];

                u8 r = static_cast<u8>(255 * clamp01(pixel.r_));
                u8 g = static_cast<u8>(255 * clamp01(pixel.g_));
                u8 b = static_cast<u8>(255 * clamp01(pixel.b_));
                fprintf(file, "%d %d %d ", r, g, b);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }

    Vector3 normalize(const Vector3& v)
    {
        f32 il = 1.0f/v.length();
        return Vector3(v.x_*il, v.y_*il, v.z_*il);
    }

    f32 dot(const Vector3& v0, const Vector3& v1)
    {
        return v0.x_*v1.x_ + v0.y_*v1.y_ + v0.z_*v1.z_;
    }

    Vector3 cross(const Vector3& v0, const Vector3& v1)
    {
        f32 x = v0.y_*v1.z_ - v0.z_*v1.y_;
        f32 y = v0.z_*v1.x_ - v0.x_*v1.z_;
        f32 z = v0.x_*v1.y_ - v0.y_*v1.x_;
        return Vector3(x, y, z);
    }

namespace
{
    s32 separateBy2(s32 x)
    {
        x = (x | (x << 8) | (x<<16)) & 0x0300F00FU;
        x = (x | (x << 4) | (x<< 8)) & 0x030C30C3U;
        x = (x | (x << 2) | (x<< 4)) & 0x09249249U;
        return x;
    }

    s32 combineBy2(s32 x)
    {
        x &= 0x09249249U;
        x = (x | (x>>2) | (x>> 4)) & 0x030C30C3U;
        x = (x | (x>>4) | (x>> 8)) & 0x0300F00FU;
        x = (x | (x>>8) | (x>>16)) & 0x3FFU;
        return x;
    }
}
    s32 mortonCode3(s32 x, s32 y, s32 z)
    {
        return separateBy2(x) | (separateBy2(y) << 1) | (separateBy2(z) << 2);
    }

    void rmortonCode3(s32& x, s32& y, s32& z, s32 w)
    {
        x = combineBy2(w);
        y = combineBy2(w>>1);
        z = combineBy2(w>>2);
    }


    void AABB::extend(const AABB& bbox)
    {
        bmin_.x_ = minimum(bmin_.x_, bbox.bmin_.x_);
        bmin_.y_ = minimum(bmin_.y_, bbox.bmin_.y_);
        bmin_.z_ = minimum(bmin_.z_, bbox.bmin_.z_);

        bmax_.x_ = maximum(bmax_.x_, bbox.bmax_.x_);
        bmax_.y_ = maximum(bmax_.y_, bbox.bmax_.y_);
        bmax_.z_ = maximum(bmax_.z_, bbox.bmax_.z_);
    }

    s32 AABB::maxExtentAxis() const
    {
        Vector3 extent = bmax_ - bmin_;
        s32 axis = (extent.x_<extent.y_)? 1 : 0;
        axis = (extent.z_<extent[axis])? axis : 2;
        return axis;
    }

    f32 AABB::halfArea() const
    {
        f32 dx = bmax_.x_ - bmin_.x_;
        f32 dy = bmax_.y_ - bmin_.y_;
        f32 dz = bmax_.z_ - bmin_.z_;
        return (dx*dy + dy*dz + dz*dx);
    }

    bool AABB::testRay(f32& tmin, f32& tmax, const Ray& ray) const
    {
        tmin = 0.0f;
        tmax = ray.t_;

        for(s32 i=0; i<3; ++i){
            if(absolute(ray.direction_[i])<F32_HITEPSILON){
                //線分とスラブが平行で、原点がスラブの中にない
                if(ray.origin_[i]<bmin_[i] || bmax_[i]<ray.origin_[i]){
                    return false;
                }

            }else{
                f32 invD = ray.invDirection_[i];
                f32 t1 = (bmin_[i] - ray.origin_[i]) * invD;
                f32 t2 = (bmax_[i] - ray.origin_[i]) * invD;

                if(t1>t2){
                    if(t2>tmin) tmin = t2;
                    if(t1<tmax) tmax = t1;
                }else{
                    if(t1>tmin) tmin = t1;
                    if(t2<tmax) tmax = t2;
                }

                if(tmin > tmax){
                    return false;
                }
                if(tmax < 0.0f){
                    return false;
                }
            }
        }
        return true;
    }

    f32 Point::getCentroidX() const
    {
        return position_.x_;
    }

    f32 Point::getCentroidY() const
    {
        return position_.y_;
    }

    f32 Point::getCentroidZ() const
    {
        return position_.z_;
    }

    Vector3 Point::getCentroid() const
    {
        return position_;
    }

    f32 Point::getCentroid(s32 axis) const
    {
        return position_[axis];
    }

    AABB Point::getBBox() const
    {
        return AABB(position_, position_);
    }

    //--- Face
    //--------------------------------------------------------------
    f32 Face::getCentroidX() const
    {
        f32 bmin = point_[0].x_;
        f32 bmax = point_[0].x_;
        for(s32 i=1; i<3; ++i){
            bmin = minimum(bmin, point_[i].x_);
            bmax = maximum(bmax, point_[i].x_);
        }
        return 0.5f * (bmin + bmax);
    }

    f32 Face::getCentroidY() const
    {
        f32 bmin = point_[0].y_;
        f32 bmax = point_[0].y_;
        for(s32 i=1; i<3; ++i){
            bmin = minimum(bmin, point_[i].y_);
            bmax = maximum(bmax, point_[i].y_);
        }
        return 0.5f * (bmin + bmax);
    }

    f32 Face::getCentroidZ() const
    {
        f32 bmin = point_[0].z_;
        f32 bmax = point_[0].z_;
        for(s32 i=1; i<3; ++i){
            bmin = minimum(bmin, point_[i].z_);
            bmax = maximum(bmax, point_[i].z_);
        }
        return 0.5f * (bmin + bmax);
    }

    Vector3 Face::getCentroid() const
    {
        Vector3 p = point_[0] + point_[1] + point_[2];
        return p*(1.0f/3.0f);
    }

    f32 Face::getCentroid(s32 axis) const
    {
        f32 p = point_[0][axis] + point_[1][axis] + point_[2][axis];
        return p*(1.0f/3.0f);
    }

    AABB Face::getBBox() const
    {
        Vector3 bmin = point_[0];
        Vector3 bmax = point_[0];
        for(s32 i=1; i<3; ++i){
            bmin.x_ = minimum(bmin.x_, point_[i].x_);
            bmin.y_ = minimum(bmin.y_, point_[i].y_);
            bmin.z_ = minimum(bmin.z_, point_[i].z_);

            bmax.x_ = maximum(bmax.x_, point_[i].x_);
            bmax.y_ = maximum(bmax.y_, point_[i].y_);
            bmax.z_ = maximum(bmax.z_, point_[i].z_);
        }
        return AABB(bmin, bmax);
    }

    bool Face::testRay(f32& t, const Ray& ray) const
    {
        Vector3 d0 = point_[1]-point_[0];
        Vector3 d1 = point_[2]-point_[0];
        Vector3 c = cross(ray.direction_, d1);

        Vector3 tvec;
        f32 discr = dot(c, d0);
        Vector3 qvec;
        if(F32_EPSILON<discr){
            //表面
            tvec = ray.origin_-point_[0];
            f32 v = dot(tvec, c);
            if(v<0.0f || discr<v){
                return false;
            }
            qvec = cross(tvec, d0);
            f32 w = dot(qvec, ray.direction_);
            if(w<0.0f || discr<(v+w)){
                return false;
            }


        } else if(discr < -F32_EPSILON){
            //裏面
            tvec = ray.origin_-point_[0];
            f32 v = dot(tvec, c);
            if(0.0f<v || v<discr){
                return false;
            }
            qvec = cross(tvec, d0);
            f32 w = dot(qvec, ray.direction_);
            if(0.0f<w || (v+w)<discr){
                return false;
            }

        } else{
            return false;
        }

        f32 invDiscr = 1.0f/discr;

        t = dot(d1, qvec);
        t *= invDiscr;
        return true;
    }

    Vector3 Face::getNormal() const
    {
        Vector3 d0 = point_[1]-point_[0];
        Vector3 d1 = point_[2]-point_[0];
        Vector3 n = cross(d0, d1);
        return normalize(n);
    }

    //--- Sphere
    //--------------------------------------------------------------
    f32 Sphere::getCentroidX() const
    {
        return position_.x_;
    }

    f32 Sphere::getCentroidY() const
    {
        return position_.y_;
    }

    f32 Sphere::getCentroidZ() const
    {
        return position_.z_;
    }

    Vector3 Sphere::getCentroid() const
    {
        return position_;
    }

    f32 Sphere::getCentroid(s32 axis) const
    {
        return position_[axis];
    }

    AABB Sphere::getBBox() const
    {
        Vector3 bmin(position_.x_-radius_, position_.y_-radius_, position_.z_-radius_);
        Vector3 bmax(position_.x_+radius_, position_.y_+radius_, position_.z_+radius_);
        return AABB(bmin, bmax);
    }

    bool Sphere::testRay(f32& t, const Ray& ray) const
    {
        Vector3 m;
        m.x_ = ray.origin_.x_ - position_.x_;
        m.y_ = ray.origin_.y_ - position_.y_;
        m.z_ = ray.origin_.z_ - position_.z_;


        f32 b = dot(m, ray.direction_);
        f32 c = dot(m, m) - radius_ * radius_;

        // 線分の起点が球の外で、向きが球の方向と逆
        if(0.0f<c && 0.0f<b){
            return false;
        }

        f32 discr = b*b - c; //判別式
        if(discr < 0.0f){
            return false;
        }

        discr = ::sqrtf(discr);
        b = -b;

        t = b - discr;
        f32 tmax = b + discr;
        return (c<=0.0f)? true : (tmax<=ray.t_);
    }

    Ray::Ray(const Vector3& origin,
        const Vector3& direction,
        f32 t)
        :origin_(origin)
        ,direction_(direction)
        ,t_(t)
    {
        invertDirection();
    }

    void Ray::invertDirection()
    {
        for(s32 i=0; i<3; ++i){
            if(0.0f<=direction_[i]){
                invDirection_[i] = (isZeroPositive(direction_[i]))? FLT_MAX : 1.0f/direction_[i];
            }else{
                invDirection_[i] = (isZeroNegative(direction_[i]))? -FLT_MAX : 1.0f/direction_[i];
            }
        }
    }

    void Ray::setDirection(const Vector3& direction)
    {
        direction_ = direction;
        invertDirection();
    }

    void Ray::setDirection(const Vector3& direction, const Vector3& invDirection)
    {
        direction_ = direction;
        invDirection_ = invDirection;
    }

    //---------------------------------------------------------
    //---
    //--- タイム関係
    //---
    //---------------------------------------------------------
    void sleep(u32 milliSeconds)
    {
#if defined(_WIN32) || defined(_WIN64)
        ::Sleep(milliSeconds);
#else
        timespec ts;
        ts.tv_sec = 0;
        while(1000<milliSeconds){
            ts.tv_sec += 1;
            milliSeconds -= 1000;
        }
        ts.tv_nsec = 1000000L * milliSeconds;
        nanosleep(&ts, NULL);
#endif
    }

    // カウント取得
    ClockType getPerformanceCounter()
    {
#if defined(_WIN32) || defined(_WIN64)
        LARGE_INTEGER count;
        QueryPerformanceCounter(&count);
        return count.QuadPart;
#else
        clock_t t = 0;
        t = clock();
        return t;
#endif
    }

    // 秒間カウント数
    ClockType getPerformanceFrequency()
    {
#if defined(_WIN32) || defined(_WIN64)
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        return freq.QuadPart;
#else
        return CLOCKS_PER_SEC;
#endif
    }

    // 秒単位の時間差分計算
    f64 calcTime64(ClockType prevTime, ClockType currentTime)
    {
        ClockType d = (currentTime>=prevTime)? currentTime - prevTime : std::numeric_limits<ClockType>::max() - prevTime + currentTime;
        f64 delta = static_cast<f64>(d)/getPerformanceFrequency();
        return delta;
    }


#define ACC_RANDOM_XORSHIFT_PROC \
    u32 t = x_^(x_<<11);\
    x_ = y_;\
    y_ = z_;\
    z_ = w_;\
    w_ = (w_^(w_>>19)) ^ (t^(t>>8));\

    RandomXorshift::RandomXorshift()
        :x_(123459876)
        ,y_(362436069)
        ,z_(521288629)
        ,w_(88675123)
    {
    }

    RandomXorshift::RandomXorshift(u32 seed)
    {
        srand(seed);
    }

    RandomXorshift::~RandomXorshift()
    {
    }

    void RandomXorshift::srand(u32 seed)
    {
        x_ = seed & 0xFFFFFFFFU;
        y_ = rand(x_, 1);
        z_ = rand(y_, 2);
        w_ = rand(z_, 3);
    }

    u32 RandomXorshift::rand()
    {
ACC_RANDOM_XORSHIFT_PROC
        return w_;
    }

    f32 RandomXorshift::frand()
    {
ACC_RANDOM_XORSHIFT_PROC

        static const u32 m0 = 0x3F800000U;
        static const u32 m1 = 0x007FFFFFU;
        t = m0|(w_&m1);
        return (*(f32*)&t) - 0.999999881f;
    }

    f32 RandomXorshift::frand2()
    {
ACC_RANDOM_XORSHIFT_PROC

        static const u32 m0 = 0x3F800000U;
        static const u32 m1 = 0x007FFFFFU;
        t = m0|(w_&m1);
        return (*(f32*)&t) - 1.000000000f;
    }

    f64 RandomXorshift::drand()
    {
        return rand()*(1.0/4294967295.0); 
    }

    u32 RandomXorshift::rand(u32 v, u32 i)
    {
        return (1812433253 * (v^(v >> 30)) + i);
    }
#undef ACC_RANDOM_XORSHIFT_PROC


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
        const __m128 bbox[2][3])
    {
        for(s32 i=0; i<3; ++i){
            tmin = _mm_max_ps(
                tmin,
                _mm_mul_ps(_mm_sub_ps(bbox[sign[i]][i], origin[i]), invDir[i]));

            tmax = _mm_min_ps(
                tmax,
                _mm_mul_ps(_mm_sub_ps(bbox[1-sign[i]][i], origin[i]), invDir[i]));
        }

        return _mm_movemask_ps(_mm_cmpge_ps(tmax, tmin));
    }

    //AABBの交差判定
    s32 testAABB(const __m128 bbox0[2][3], const __m128 bbox1[2][3])
    {
        u32 mask = 0xFFFFFFFFU;
        f32 fmask = *((f32*)&mask);

        __m128 t = _mm_set1_ps(fmask);
        for(s32 i=0; i<3; ++i){
            t = _mm_and_ps(t, _mm_cmple_ps(bbox0[0][i], bbox1[1][i]));
            t = _mm_and_ps(t, _mm_cmple_ps(bbox1[0][i], bbox0[1][i]));
        }
        return _mm_movemask_ps(t);
    }

    //
    s32 testSphereAABB(const __m128 position[3], const __m128 radius, const __m128 bbox[2][3])
    {
        u32 mask = 0xFFFFFFFFU;
        f32 fmask = *((f32*)&mask);

        __m128 tbbox[2][3];
        for(s32 i=0; i<3; ++i){
            tbbox[0][i] = _mm_sub_ps(bbox[0][i], radius);
            tbbox[1][i] = _mm_add_ps(bbox[1][i], radius);
        }
        __m128 t = _mm_set1_ps(fmask);
        for(s32 i=0; i<3; ++i){
            t = _mm_and_ps(t, _mm_cmple_ps(position[i], tbbox[1][i]));
            t = _mm_and_ps(t, _mm_cmple_ps(tbbox[0][i], position[i]));
        }
        return _mm_movemask_ps(t);
    }
}
}
