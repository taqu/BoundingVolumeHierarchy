/**
@file lrender.cpp
@author t-sakai
@date 2015/09/04 create
*/
#include "lrender.h"
#include <stdio.h>

namespace lrender
{
    const f32 Epsilon = 1.0e-5f;
    const f32 HitEpsilon = 1.0e-5f;

    void printImage(const char* filename, RGB* rgb, s32 width, s32 height)
    {
        LASSERT(NULL != filename);
        LASSERT(NULL != rgb);

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

                u8 r = 255 * clamp01(pixel.r_);
                u8 g = 255 * clamp01(pixel.g_);
                u8 b = 255 * clamp01(pixel.b_);
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

    void BBox::extend(const BBox& bbox)
    {
        bmin_.x_ = minimum(bmin_.x_, bbox.bmin_.x_);
        bmin_.y_ = minimum(bmin_.y_, bbox.bmin_.y_);
        bmin_.z_ = minimum(bmin_.z_, bbox.bmin_.z_);

        bmax_.x_ = maximum(bmax_.x_, bbox.bmax_.x_);
        bmax_.y_ = maximum(bmax_.y_, bbox.bmax_.y_);
        bmax_.z_ = maximum(bmax_.z_, bbox.bmax_.z_);
    }

    s32 BBox::maxExtentAxis() const
    {
        Vector3 extent = bmax_ - bmin_;
        s32 axis = (extent.x_<extent.y_)? 1 : 0;
        axis = (extent.z_<extent[axis])? axis : 2;
        return axis;
    }

    f32 BBox::halfArea() const
    {
        f32 dx = bmax_.x_ - bmin_.x_;
        f32 dy = bmax_.y_ - bmin_.y_;
        f32 dz = bmax_.z_ - bmin_.z_;
        return (dx*dy + dy*dz + dz*dx);
    }

    bool BBox::testRay(f32& tmin, f32& tmax, const Ray& ray) const
    {
        tmin = 0.0f;
        tmax = ray.t_;

        for(s32 i=0; i<3; ++i){
            if(absolute(ray.direction_[i])<HitEpsilon){
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

    BBox Point::getBBox() const
    {
        return BBox(position_, position_);
    }


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

    BBox Face::getBBox() const
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
        return BBox(bmin, bmax);
    }

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

    BBox Sphere::getBBox() const
    {
        Vector3 bmin(position_.x_-radius_, position_.y_-radius_, position_.z_-radius_);
        Vector3 bmax(position_.x_+radius_, position_.y_+radius_, position_.z_+radius_);
        return BBox(bmin, bmax);
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

    //void HitRecord::calcPoint(Vector3& dst)
    //{
    //    Vector3 v0 = (*face_)[0];
    //    Vector3 v1 = (*face_)[1];
    //    Vector3 v2 = (*face_)[2];

    //    v0 *= (1.0f - u_ - v_);
    //    v1 *= u_;
    //    v2 *= v_;

    //    dst = v0;
    //    dst += v1;
    //    dst += v2;
    //}

    //void HitRecord::calcNormal(Vector3& normal)
    //{
    //    Vector3 v0 = face_->getNormal(0);
    //    Vector3 v1 = face_->getNormal(1);
    //    Vector3 v2 = face_->getNormal(2);

    //    v0 *= (1.0f - u_ - v_);
    //    v1 *= u_;
    //    v2 *= v_;

    //    normal = v0;
    //    normal += v1;
    //    normal += v2;
    //    normal.normalize();
    //}

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

    // ミリ秒単位の時間を取得
    u32 getTime()
    {
#if defined(_WIN32) || defined(_WIN64)
        DWORD time = timeGetTime();
        return static_cast<u32>(time);
#else
        struct timeval tv;
        gettimeofday(&tv, NULL);

        return static_cast<u32>(tv.tv_sec*1000 + tv.tv_usec/1000);
#endif
    }

#define LCORE_RANDOM_XORSHIFT_PROC \
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
LCORE_RANDOM_XORSHIFT_PROC
        return w_;
    }

    f32 RandomXorshift::frand()
    {
LCORE_RANDOM_XORSHIFT_PROC

        static const u32 m0 = 0x3F800000U;
        static const u32 m1 = 0x007FFFFFU;
        t = m0|(w_&m1);
        return (*(f32*)&t) - 0.999999881f;
    }

    f32 RandomXorshift::frand2()
    {
LCORE_RANDOM_XORSHIFT_PROC

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
#undef LCORE_RANDOM_XORSHIFT_PROC


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
