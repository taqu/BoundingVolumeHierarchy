#include "BinBVH.h"
#include "MedBVH.h"
#include "MedQBVH.h"
#include "BinQBVH.h"
#include "LQBVH.h"
#include "HLQBVH.h"
#include "GRIDQBVH.h"
#include "sibenik/sibenik.bin"
#include "conference/conference.bin"

using namespace accel;

namespace
{
    s32 width = 640;
    s32 height = 480;

    Vector3 cameraPosition;
    Vector3 cameraForward;
    Vector3 cameraUp;
    Vector3 cameraRight;
    f32 iw;
    f32 ih;

    Vector3 sx;
    Vector3 sy;

    RGB* rgb = NULL;

    void initScene(const Vector3 position, const Vector3& direction)
    {
        cameraPosition = position;
        cameraForward = normalize(direction);
        cameraUp = Vector3(0.0f, 1.0f, 0.0f);
        cameraRight = normalize(cross(cameraUp, cameraForward));
        f32 aspect = static_cast<f32>(width)/height;

        f32 fovy = ::tanf(0.5f * 60.0f/180.0f*static_cast<f32>(M_PI));
        f32 dx = fovy * aspect;
        f32 dy = fovy;

        iw = 1.0f/(width-1);
        ih = 1.0f/(height-1);
        sx = dx * cameraRight;
        sy = dy * cameraUp;

        rgb = (RGB*)ACC_MALLOC(sizeof(RGB)*width*height);
    }

    void termScene()
    {
        ACC_FREE(rgb);
    }

#if 0
    s32 divide(u8 axis, f32 mid, s32 start, s32 size, Face* tries)
    {
        typedef PrimitivePolicy<Face> PrimitivePolicy;

        if(size<=0){
            return 0;
        }
        s32 left = start;
        s32 right = start + size-1;
        for(;;){
            while(left<right){
                f32 l = PrimitivePolicy::getCentroid(tries[left], axis);
                if(mid<l){
                    break;
                }
                ++left;
            }
            while(left<right){
                f32 r = PrimitivePolicy::getCentroid(tries[right], axis);
                if(r<=mid){
                    break;
                }
                --right;
            }
            if(right<=left){
                break;
            }
            swap(tries[left], tries[right]);
            ++left;
            --right;
        }
        return right;
    }
#endif
}

void brute_force(
    s32& depth,
    f64& buildTime,
    f64& renderTime,
    f64& raysPerSecond,
    s32 numPrimitives,
    Face* faces,
    const char* filename)
{
    Timer<true> timer;

    timer.start();
    buildTime = 0.0f;
    timer.stop();
    buildTime = timer.getAverage();

    depth = 0;
    timer.reset();

    Ray ray;
    ray.origin_ = cameraPosition;
    ray.t_ = 1.0e32f;

    timer.start();
    for(s32 i=0; i<height; ++i){
        Vector3 vy = sy * (1.0f - 2.0f*i*ih);
        for(s32 j=0; j<width; ++j){
            Vector3 vx = sx * (2.0f*j*iw - 1.0f);
            
            ray.direction_ = normalize(cameraForward+vx+vy);
            ray.invertDirection();
            ray.t_ = 1.0e32f;
            RGB& pixel = rgb[i*width+j];

            Face* face = NULL;
            for(s32 k=0; k<numPrimitives; ++k){
                f32 t;
                if(faces[k].testRay(t, ray)){
                    if(t<ray.t_){
                        ray.t_ = t;
                        face = &faces[k];
                    }
                }
            }
            if(NULL != face){
                Vector3 normal = face->getNormal();
                pixel.r_ = pixel.g_ = pixel.b_ = maximum(normal.y_, 0.0f);
            }else{
                pixel.r_ = pixel.g_ = pixel.b_ = 0;
            }
        }
    }
    timer.stop();
    renderTime = timer.getAverage();
    raysPerSecond = (height*width)/renderTime;
    printImage(filename, rgb, width, height);
    printf("%s is done.\n", filename);
}

void brute_force(
    FILE* file,
    s32 numPrimitives,
    Face* faces,
    const char* name,
    const char* filename)
{
    s32 depth=0;
    f64 buildTime=0.0;
    f64 renderTime=0.0;
    f64 raysPerSecond=0.0;
    brute_force(depth, buildTime, renderTime, raysPerSecond, numPrimitives, faces, filename);
    fprintf(file, "%s depth: %d, build: %f, render: %f, rays: %f\r\n", name, depth, buildTime, renderTime, raysPerSecond);
}

template<class T>
void test(
    s32& depth,
    f64& buildTime,
    f64& renderTime,
    f64& raysPerSecond,
    s32 numPrimitives,
    Face* faces,
    const char* filename)
{
    Timer<true> timer;
    T* bvh = new T;

    timer.start();
    bvh->build(numPrimitives, faces);
    timer.stop();
    buildTime = timer.getAverage();

    depth = bvh->getDepth();
    timer.reset();

    Ray ray;
    ray.origin_ = cameraPosition;
    ray.t_ = 1.0e32f;

    timer.start();
    for(s32 i=0; i<height; ++i){
        Vector3 vy = sy * (1.0f - 2.0f*i*ih);
        for(s32 j=0; j<width; ++j){
            Vector3 vx = sx * (2.0f*j*iw - 1.0f);

            ray.direction_ = normalize(cameraForward+vx+vy);
            ray.invertDirection();
            ray.t_ = 1.0e32f;
            RGB& pixel = rgb[i*width+j];
            HitRecord hitRecord = bvh->intersect(ray);
            if(NULL != hitRecord.primitive_){
                Vector3 normal = reinterpret_cast<const Face*>(hitRecord.primitive_)->getNormal();
                normal.y_ = absolute(normal.y_);
                pixel.r_ = pixel.g_ = pixel.b_ = maximum(normal.y_, 0.0f);
            }else{
                pixel.r_ = pixel.g_ = pixel.b_ = 0;
            }
        }
    }
    timer.stop();
    delete bvh;
    renderTime = timer.getAverage();
    raysPerSecond = (height*width)/renderTime;
    printImage(filename, rgb, width, height);
    printf("%s is done.\n", filename);
}

template<class T>
void test(
    FILE* file,
    s32 numPrimitives,
    Face* faces,
    const char* name,
    const char* filename)
{
    s32 depth=0;
    f64 buildTime=0.0;
    f64 renderTime=0.0;
    f64 raysPerSecond=0.0;
    test<T>(depth, buildTime, renderTime, raysPerSecond, numPrimitives, faces, filename);
    fprintf(file, "%s depth: %d, build: %f, render: %f, rays: %f\r\n", name, depth, buildTime, renderTime, raysPerSecond);
}

int main(int /*argc*/, char** /*argv*/)
{
    FILE* file = NULL;
    fopen_s(&file, "statistics.txt", "wb");
    if(NULL == file){
        return 0;
    }

#if 1
//#define ENABLE_BRUTE_FORCE (1)
    const int SibenikNumFaces = sizeof(sibenik_indices)/sizeof(sibenik_indices[0])/3;
    Face* faces = (Face*)ACC_MALLOC(sizeof(Face) * SibenikNumFaces);
    for(s32 i=0; i<SibenikNumFaces; ++i){
        s32 index = i*3;
        s32 idx0 = sibenik_indices[index+0]*3;
        s32 idx1 = sibenik_indices[index+1]*3;
        s32 idx2 = sibenik_indices[index+2]*3;

        faces[i].point_[0] = Vector3(sibenik_vertices[idx0+0], sibenik_vertices[idx0+1], sibenik_vertices[idx0+2]);
        faces[i].point_[1] = Vector3(sibenik_vertices[idx1+0], sibenik_vertices[idx1+1], sibenik_vertices[idx1+2]);
        faces[i].point_[2] = Vector3(sibenik_vertices[idx2+0], sibenik_vertices[idx2+1], sibenik_vertices[idx2+2]);
    }
    initScene(Vector3(-15.0f, -5.0f, 0.0f), Vector3(1.0f, -0.2f, 0.0f));

    fprintf(file, "sibenik num primitives: %d\r\n", SibenikNumFaces);
#ifdef ENABLE_BRUTE_FORCE
    s32 brute_force_depth = 0;
    f64 brute_force_build_time=0.0, brute_force_RaysPerSecond=0.0;
    brute_force(file, SibenikNumFaces, faces, "brute force", "sibenik_brute_force.ppm");
#endif
    test<MedBVH<Face> >(file, SibenikNumFaces, faces, "medbvh", "sibenik_medbvh.ppm");

    test<BinBVH<Face> >(file, SibenikNumFaces, faces, "binbvh", "sibenik_binbvh.ppm");

    test<MedQBVH<Face> >(file, SibenikNumFaces, faces, "medqbvh", "sibenik_medqbvh.ppm");

    test<BinQBVH<Face> >(file, SibenikNumFaces, faces, "binqbvh", "sibenik_binqbvh.ppm");

    //test<LQBVH<Face> >(file, SibenikNumFaces, faces, "lqbvh", "sibenik_lqbvh.ppm");

    //test<HLQBVH<Face> >(file, SibenikNumFaces, faces, "hlqbvh", "sibenik_hlqbvh.ppm");

    test<GRIDQBVH<Face> >(file, SibenikNumFaces, faces, "gridqbvh", "sibenik_gridqbvh.ppm");

    termScene();

    ACC_FREE(faces);


    const int ConferenceNumFaces = sizeof(conference_indices)/sizeof(conference_indices[0])/3;
    faces = (Face*)ACC_MALLOC(sizeof(Face) * ConferenceNumFaces);
    for(s32 i=0; i<ConferenceNumFaces; ++i){
        s32 index = i*3;
        s32 idx0 = conference_indices[index+0]*3;
        s32 idx1 = conference_indices[index+1]*3;
        s32 idx2 = conference_indices[index+2]*3;

        faces[i].point_[0] = Vector3(conference_vertices[idx0+0], conference_vertices[idx0+1], conference_vertices[idx0+2]);
        faces[i].point_[1] = Vector3(conference_vertices[idx1+0], conference_vertices[idx1+1], conference_vertices[idx1+2]);
        faces[i].point_[2] = Vector3(conference_vertices[idx2+0], conference_vertices[idx2+1], conference_vertices[idx2+2]);
    }
    initScene(Vector3(1610.0f, 356.0f, -322.0f), Vector3(-100.0f, 0.0f, 0.0f));

    fprintf(file, "conference num primitives: %d\r\n", ConferenceNumFaces);
#ifdef ENABLE_BRUTE_FORCE
    s32 brute_force_depth = 0;
    f64 brute_force_build_time=0.0, brute_force_RaysPerSecond=0.0;
    brute_force(file, ConferenceNumFaces, faces, "brute force", "conference_brute_force.ppm");
#endif
    test<MedBVH<Face> >(file, ConferenceNumFaces, faces, "medbvh", "conference_medbvh.ppm");

    test<BinBVH<Face> >(file, ConferenceNumFaces, faces, "binbvh", "conference_binbvh.ppm");

    test<MedQBVH<Face> >(file, ConferenceNumFaces, faces, "medqbvh", "conference_medqbvh.ppm");

    test<BinQBVH<Face> >(file, ConferenceNumFaces, faces, "binqbvh", "conference_binqbvh.ppm");

    //test<LQBVH<Face> >(file, ConferenceNumFaces, faces, "lqbvh", "conference_lqbvh.ppm");

    //test<HLQBVH<Face> >(file, ConferenceNumFaces, faces, "hlqbvh", "conference_hlqbvh.ppm");

    test<GRIDQBVH<Face> >(file, ConferenceNumFaces, faces, "gridqbvh", "conference_gridqbvh.ppm");

    termScene();

    ACC_FREE(faces);

#else
//#define ENABLE_BRUTE_FORCE (1)
    RandomXorshift random;
    random.srand(123456789);

    static const f32 Area = 6.0f;
    static const int NumFaces = 32;
    Face* faces = (Face*)ACC_MALLOC(sizeof(Face) * NumFaces);
    for(s32 i=0; i<NumFaces; ++i){
        f32 px = (random.frand()*2.0f - 1.0f)*Area;
        f32 py = (random.frand()*2.0f - 1.0f)*Area;
        f32 pz = (random.frand()*2.0f - 1.0f)*Area;
        for(s32 j=0; j<3; ++j){
            faces[i].point_[j].x_ = (random.frand()*2.0f - 1.0f)*0.2f + px;
            faces[i].point_[j].y_ = (random.frand()*2.0f - 1.0f)*0.2f + py;
            faces[i].point_[j].z_ = (random.frand()*2.0f - 1.0f)*0.2f + pz;
        }
    }
    initScene(Vector3(0.0f, 0.0f, 12.0f), Vector3(0.0f, 0.0f, -1.0f));

    fprintf(file, "num primitives: %d\r\n", NumFaces);

#ifdef ENABLE_BRUTE_FORCE
    s32 brute_force_depth = 0;
    f64 brute_force_build_time=0.0, brute_force_RaysPerSecond=0.0;
    brute_force(file, NumFaces, faces, "brute force", "out_brute_force.ppm");
#endif
    test<MedBVH<Face> >(file, NumFaces, faces, "medbvh", "out_medbvh.ppm");

    test<BinBVH<Face> >(file, NumFaces, faces, "binbvh", "out_binbvh.ppm");

    test<MedQBVH<Face> >(file, NumFaces, faces, "medqbvh", "out_medqbvh.ppm");

    test<BinQBVH<Face> >(file, NumFaces, faces, "binqbvh", "out_binqbvh.ppm");

    test<LQBVH<Face> >(file, NumFaces, faces, "lqbvh", "out_lqbvh.ppm");

    //test<HLQBVH<Face> >(file, NumFaces, faces, "hlqbvh", "out_hlqbvh.ppm");

    test<GRIDQBVH<Face> >(file, NumFaces, faces, "gridqbvh", "out_gridqbvh.ppm");

    termScene();

    ACC_FREE(faces);

#endif

    fclose(file);
    return 0;
}
