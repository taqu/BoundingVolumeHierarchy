#include "BVH.h"
#include "BinBVH.h"
#include "MidBVH.h"
#include "MidQBVH.h"
#include "BinQBVH.h"

namespace
{
    lrender::s32 width = 640;
    lrender::s32 height = 480;

    lrender::Vector3 cameraPosition;
    lrender::Vector3 cameraForward;
    lrender::Vector3 cameraUp;
    lrender::Vector3 cameraRight;

    lrender::f32 iw;
    lrender::f32 ih;

    lrender::Vector3 sx;
    lrender::Vector3 sy;

    lrender::RGB* rgb = NULL;

    void initScene()
    {
        cameraPosition = lrender::Vector3(0.0f, 0.0f, 20.0f);
        cameraForward = lrender::Vector3(0.0f, 0.0f, -1.0f);
        cameraUp = lrender::Vector3(0.0f, 1.0f, 0.0f);
        cameraRight = lrender::normalize(lrender::cross(cameraUp, cameraForward));
        lrender::f32 aspect = static_cast<lrender::f32>(width)/height;

        lrender::f32 fovy = ::tanf(0.5f * 60.0f/180.0f*static_cast<lrender::f32>(M_PI));
        lrender::f32 dx = fovy * aspect;
        lrender::f32 dy = fovy;

        iw = 2.0f/(width-1);
        ih = 2.0f/(height-1);
        sx = dx * cameraRight;
        sy = dy * cameraUp;

        rgb = (lrender::RGB*)LIME_MALLOC(sizeof(lrender::RGB)*width*height);
    }

    void termScene()
    {
        LIME_FREE(rgb);
    }
}

template<class T>
void test(
    lrender::s32& depth,
    lrender::f64& buildTime,
    lrender::f64& renderTime,
    lrender::s32 numPrimitives,
    lrender::Sphere* spheres,
    const char* filename,
    bool debugOut)
{
    lrender::Timer<true> timer;
    T bvh;

    timer.begin();
    bvh.build(numPrimitives, spheres);
    timer.end();
    buildTime = timer.getAverage();

    if(debugOut){
        bvh.print("bvh.txt");
    }

    depth = bvh.getDepth();
    timer.reset();

    lrender::Ray ray;
    ray.origin_ = cameraPosition;
    ray.t_ = 1.0e32f;

    timer.begin();
    for(lrender::s32 i=0; i<height; ++i){
        lrender::Vector3 vy = sy * (1.0f - i*ih);
        for(lrender::s32 j=0; j<width; ++j){
            lrender::Vector3 vx = sx * (j*iw - 1.0f);
            
            ray.direction_ = lrender::normalize(cameraForward+vx+vy);
            ray.invertDirection();
            ray.t_ = 1.0e32f;
            lrender::RGB& pixel = rgb[i*width+j];
            lrender::HitRecord hitRecord = bvh.intersect(ray);
            if(NULL != hitRecord.primitive_){
                lrender::u8 d = static_cast<lrender::u8>(255.0f*hitRecord.t_/200.0f);
                pixel.r_ = pixel.g_ = pixel.b_ = d;
            }else{
                pixel.r_ = pixel.g_ = pixel.b_ = 0;
            }
        }
    }
    timer.end();

    renderTime = timer.getAverage();
    lrender::printImage(filename, rgb, width, height);
}

int main(int argc, char** argv)
{
    lrender::RandomXorshift random;
    random.srand(123456789);

    static const lrender::f32 Radius = 0.2f;
    static const lrender::f32 Area = 10.0f;

#if 0
    static const int NumSpheres = 8;
    lrender::Sphere* spheres = (lrender::Sphere*)LIME_MALLOC(sizeof(lrender::Sphere) * NumSpheres);

    spheres[0] = lrender::Sphere(lrender::Vector3(-Area,  Area, Area), Radius);
    spheres[1] = lrender::Sphere(lrender::Vector3( Area,  Area, Area), Radius);
    spheres[2] = lrender::Sphere(lrender::Vector3(-Area, -Area, Area), Radius);
    spheres[3] = lrender::Sphere(lrender::Vector3( Area, -Area, Area), Radius);

    spheres[4] = lrender::Sphere(lrender::Vector3(-Area,  Area, -Area), Radius);
    spheres[5] = lrender::Sphere(lrender::Vector3( Area,  Area, -Area), Radius);
    spheres[6] = lrender::Sphere(lrender::Vector3(-Area, -Area, -Area), Radius);
    spheres[7] = lrender::Sphere(lrender::Vector3( Area, -Area, -Area), Radius);
#else
    static const int NumSpheres = 100*1024;
    lrender::Sphere* spheres = (lrender::Sphere*)LIME_MALLOC(sizeof(lrender::Sphere) * NumSpheres);
    for(lrender::s32 i=0; i<NumSpheres; ++i){
        lrender::f32 px = (random.frand()*2.0f - 1.0f)*Area;
        lrender::f32 py = (random.frand()*2.0f - 1.0f)*Area;
        lrender::f32 pz = (random.frand()*2.0f - 1.0f)*Area;
        lrender::f32 r = (random.frand() + 0.5f)*Radius;
        spheres[i] = lrender::Sphere(lrender::Vector3(px, py, pz), r);
    }
#endif
    initScene();

    lrender::s32 bvhDepth;
    lrender::f64 bvhBuildTime, bvhRenderTime;
    test<lrender::BVH<lrender::Sphere> >(bvhDepth, bvhBuildTime, bvhRenderTime, NumSpheres, spheres, "out_bvh.ppm", false);

    lrender::s32 midbvhDepth;
    lrender::f64 midbvhBuildTime, midbvhRenderTime;
    test<lrender::MidBVH<lrender::Sphere> >(midbvhDepth, midbvhBuildTime, midbvhRenderTime, NumSpheres, spheres, "out_midbvh.ppm", false);

    lrender::s32 binbvhDepth;
    lrender::f64 binbvhBuildTime, binbvhRenderTime;
    test<lrender::BinBVH<lrender::Sphere> >(binbvhDepth, binbvhBuildTime, binbvhRenderTime, NumSpheres, spheres, "out_binbvh.ppm", false);

    lrender::s32 qbvhDepth;
    lrender::f64 qbvhBuildTime, qbvhRenderTime;
    test<lrender::MidQBVH<lrender::Sphere> >(qbvhDepth, qbvhBuildTime, qbvhRenderTime, NumSpheres, spheres, "out_qbvh.ppm", false);

    lrender::s32 binqbvhDepth;
    lrender::f64 binqbvhBuildTime, binqbvhRenderTime;
    test<lrender::BinQBVH<lrender::Sphere> >(binqbvhDepth, binqbvhBuildTime, binqbvhRenderTime, NumSpheres, spheres, "out_binqbvh.ppm", false);

    termScene();

    LIME_FREE(spheres);

    FILE* file = NULL;
    fopen_s(&file, "statistics.txt", "wb");
    if(NULL != file){
        fprintf(file, "num primitives: %d\r\n", NumSpheres);
        fprintf(file, "sahbvh depth: %d, build: %f, render: %f\r\n", bvhDepth, bvhBuildTime, bvhRenderTime);
        fprintf(file, "midbvh depth: %d, build: %f, render: %f\r\n", midbvhDepth, midbvhBuildTime, midbvhRenderTime);
        fprintf(file, "binbvh depth: %d, build: %f, render: %f\r\n", binbvhDepth, binbvhBuildTime, binbvhRenderTime);
        fprintf(file, "\r\n\r\n");
        fprintf(file, "midqbvh depth: %d, build: %f, render: %f\r\n", qbvhDepth, qbvhBuildTime, qbvhRenderTime);
        fprintf(file, "binqbvh depth: %d, build: %f, render: %f\r\n", binqbvhDepth, binqbvhBuildTime, binqbvhRenderTime);
        fclose(file);
    }
    return 0;
}
