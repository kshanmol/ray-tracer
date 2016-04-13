#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>

#define BLOCK_SIZE 32
#define HD __host__ __device__
#define WIDTH 1024
#define REFLECTION_DEPTH 1

static const double eps = 1e-8;

HD double det(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3);

int clamp(int what, int low, int high);

template<typename T>
class Vec3{
public:

    T x, y, z;
    //Constructors
    HD Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    HD Vec3(T val) : x(val), y(val), z(val) {}
    HD Vec3(T xval, T yval, T zval) : x(xval), y(yval), z(zval) {}

    HD Vec3& normalize(){
        T nor2 = length2();
        if (nor2 > 0) {
            T nor_inv = 1 / sqrt(nor2);
            x *= nor_inv, y *= nor_inv, z *= nor_inv;
        }
        return *this;
    }

    HD T dotProduct(const Vec3<T> &v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    HD Vec3<T> crossProduct(const Vec3<T> &v) const {

        T tmpX = y * v.z - z * v.y;
        T tmpY = z * v.x - x * v.z;
        T tmpZ = x * v.y - y * v.x;
        return Vec3<T>(tmpX, tmpY, tmpZ );
    }

    HD T length2(){
        return x * x + y * y + z * z;
    }
    HD T length(){
        return sqrt(length2());
    }

    HD Vec3<T> scale(const T &f) const {
        return Vec3<T>(x * f, y * f, z * f);
    }

    HD Vec3<T> multiply(const Vec3<T> &v) const {
        return Vec3<T>(x * v.x, y * v.y, z * v.z);
    }

    HD Vec3<T> subtract(const Vec3<T> &v) const {
        return Vec3<T>(x - v.x, y - v.y, z - v.z);
    }

    HD Vec3<T> add(const Vec3<T> &v) const {
        return Vec3<T>(x + v.x, y + v.y, z + v.z);
    }

    HD Vec3<T> negate() const {
        return Vec3<T>(-x, -y, -z);
    }

    //Helper function to format display
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v){
        os << "(" << v.x << " " << v.y << " " << v.z << ")";
        return os;
    }
};

typedef Vec3<float> Vec3f;

HD Vec3f reflect(const Vec3f &I, const Vec3f &N){
	return I.subtract(N.scale(2*I.dotProduct(N))).negate();
}

class Triangle{
public:

    Vec3f v0, v1, v2;
    Vec3f tv0, tv1, tv2; // texture coordinates of vertices
    bool reflective;

    HD
    Triangle(
        const Vec3f &v_0,
        const Vec3f &v_1,
        const Vec3f &v_2,
        const Vec3f &tv_0,
        const Vec3f &tv_1,
        const Vec3f &tv_2) :
        v0(v_0), v1(v_1), v2(v_2),
        tv0(tv_0), tv1(tv_1), tv2(tv_2)
    { /* empty */ }

    HD
    Triangle(
        const Vec3f &v_0,
        const Vec3f &v_1,
        const Vec3f &v_2,
        const Vec3f &tv_0,
        const Vec3f &tv_1,
        const Vec3f &tv_2,
        const bool reflect_) :
        v0(v_0), v1(v_1), v2(v_2),
        tv0(tv_0), tv1(tv_1), tv2(tv_2),
        reflective(reflect_)
    { /* empty */ }

    HD bool rayTriangleIntersect(const Vec3f &orig, const Vec3f &dir, float &t, float &beta, float &gamma){

        double A = det(
                    v0.x - v1.x, v0.x - v2.x, dir.x,
                    v0.y - v1.y, v0.y - v2.y, dir.y,
                    v0.z - v1.z, v0.z - v2.z, dir.z
                    );

        double t_ = det(
                    v0.x - v1.x, v0.x - v2.x, v0.x - orig.x,
                    v0.y - v1.y, v0.y - v2.y, v0.y - orig.y,
                    v0.z - v1.z, v0.z - v2.z, v0.z - orig.z
                    );
        t_ = t_ / A;

        double beta_ = det(
                    v0.x - orig.x, v0.x - v2.x, dir.x,
                    v0.y - orig.y, v0.y - v2.y, dir.y,
                    v0.z - orig.z, v0.z - v2.z, dir.z
                    );
        beta_ = beta_ / A;

        double gamma_ = det(
                    v0.x - v1.x, v0.x - orig.x, dir.x,
                    v0.y - v1.y, v0.y - orig.y, dir.y,
                    v0.z - v1.z, v0.z - orig.z, dir.z
                    );
        gamma_ = gamma_ / A;

        if (beta_ > 0 && gamma_ > 0 && beta_ + gamma_ < 1)
        {
            t = t_;
            beta = beta_;
            gamma = gamma_;
            return true;
        }
        return false;

    }

    HD Vec3f getNormal(Vec3f point) const
    {
        // from http://math.stackexchange.com/a/137551
        Vec3f p = point.subtract(v1);
        Vec3f q = v0.subtract(v2);

        return Vec3f(
                    p.y * q.z - p.z * q.y,
                    -1 * (p.x*q.z - p.z * q.x),
                    p.x*q.y - p.y*q.x
                    );
    }
};

HD Vec3f trace(Vec3f &rayorig, Vec3f &raydir, Triangle* triangle_list, int tl_size, Vec3f *lights, int n_lights, int depth=0);

HD float max_(float a, float b){ return (a < b) ? b : a; }

__global__
void trace_kernel (float* params, Triangle* triangle_list, Vec3f *light_positions, Vec3f* image){

    //Pixel coordinates
    int x = threadIdx.x + blockIdx.x*blockDim.x;
    int y = threadIdx.y + blockIdx.y*blockDim.y;

    //Unpack parameters
    float invWidth = params[0], invHeight = params[1];
    float fov = params[2], aspectratio = params[3];
    float angle = params[4];
    int tl_size = (int) params[5];
    int n_lights = (int) params[6];

    //http://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-ray-tracing/ray-tracing-practical-example

    //Set up camera view
    float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
    float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
    Vec3f raydir(xx, yy, 2);
    raydir.normalize();

    //Trace ray
    Vec3f camera_pos(0,0,7);
    image[y*WIDTH + x] = trace(camera_pos, raydir, triangle_list,tl_size, light_positions, n_lights);

}

HD Triangle *nearest_triangle(Vec3f rayorig, Vec3f raydir, Triangle *triangle_list, int tl_size, float &tnear, bool use_eps=false){
    //Ray triangle intersection
    tnear = INFINITY;
    float beta = INFINITY, // TODO: check if still needed
          gamma = INFINITY;

    Triangle* triangle_near = NULL;
    for (unsigned int i = 0; i < tl_size; ++i) {
        float t0 = INFINITY,
              beta_ = INFINITY,
              gamma_ = INFINITY;
        if (triangle_list[i].rayTriangleIntersect(rayorig, raydir, t0, beta_, gamma_)) {
            if (t0 < tnear) {
                if( (!use_eps) || (use_eps && t0 > eps) ){
                    tnear = t0;
                    triangle_near = &triangle_list[i];
                    beta = beta_;
                    gamma = gamma_;
                }
            }
        }
    }

    return triangle_near;
}

HD Vec3f trace(Vec3f &rayorig, Vec3f &raydir, Triangle* triangle_list, int tl_size, Vec3f *light_positions, int n_lights, int depth)
{

    //Ray triangle intersection
    float tnear = INFINITY;

    const Triangle* triangle_near = nearest_triangle(rayorig, raydir, triangle_list, tl_size, tnear);

    if (!triangle_near)
        return Vec3f(0);

    // Simple blinn phong shading
    Vec3f color(200.0);
    float kd = 0.3f;
    float ks = 0.5f;
    float spec_alpha = 4;

    Vec3f result_color(0.0f);
    Vec3f ambient(10.0f);

    Vec3f poi = rayorig.add( raydir.scale(tnear) );
    Vec3f n = triangle_near->getNormal(poi).normalize();

    if(triangle_near->reflective)
    {
        if (depth <= REFLECTION_DEPTH)
        {
            Vec3f reflected_ray = reflect(raydir, n).normalize();
            Vec3f reflected_shading = trace(poi, reflected_ray, triangle_list, tl_size, light_positions, n_lights, depth+1);
            return Vec3f(0, 0, 100).add(reflected_shading);
        }
    }

    for(unsigned int i = 0; i < n_lights; i++)
    {
        // TODO: extremely ineffecient color decision
        color = Vec3f( (i % 3 == 0) ? 200: 0, (i % 3 == 1) ? 200: 0, (i % 3 == 2) ? 200: 0);

        Vec3f light = light_positions[i];
        Vec3f l = poi.subtract(light).normalize();

        // checking if its in a shadow region.
        float tnear_temp;
        Vec3f shadowRay = l.negate().normalize();
        Vec3f poi_epsilon = poi.add(shadowRay.scale(eps));

        Triangle *object_between_light = nearest_triangle(poi_epsilon, shadowRay, triangle_list, tl_size, tnear_temp);
        // Triangle *object_between_light = nearest_triangle(poi, shadowRay, triangle_list, tl_size, tnear_temp, true);

        if (object_between_light != NULL){
            Vec3f diffuse = color.scale(kd * max_(float(0), n.dotProduct(l.normalize())));
            Vec3f specular = color.scale(ks * pow(max_(float(0), reflect(l,n).dotProduct(raydir.negate())), spec_alpha));
            result_color = result_color.add(diffuse).add(specular).add(ambient);
        }
    }
    return result_color.add(ambient);
}

void render(const std::vector<Triangle*> &triangle_list, const Vec3f *light_positions, int n_lights, const char *f_output){

    //Define image size, calculate camera view parameters
    int width = WIDTH, height = WIDTH;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 18);
    int tl_size = triangle_list.size();

    //Cuda events for timing data
    cudaEvent_t c_start, c_stop, t_start, t_stop, k_start, k_stop;
    cudaEventCreate(&c_start);
    cudaEventCreate(&c_stop);
    cudaEventCreate(&t_start);
    cudaEventCreate(&t_stop);
    cudaEventCreate(&k_start);
    cudaEventCreate(&k_stop);

    float h_params[7] = {invWidth, invHeight, fov, aspectratio, angle, tl_size, n_lights};

    Triangle* h_triangle_list = (Triangle*)malloc(tl_size*sizeof(Triangle));
    for(int i = 0;i<tl_size;i++)
        h_triangle_list[i] = *triangle_list[i];

    //Parallel program begins
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(width/dimBlock.x, height/dimBlock.y);

    float* d_params;
    Triangle* d_triangle_list;
    Vec3f *d_image;
    Vec3f *d_light_positions;

    //Copy parameters needed to set up camera view, triangle list
    cudaMalloc(&d_params, 7*sizeof(float));
    cudaMalloc(&d_triangle_list, tl_size*sizeof(Triangle));
    cudaMalloc(&d_image, width*height*sizeof(Vec3f));
    cudaMalloc(&d_light_positions, n_lights*sizeof(Vec3f));

    cudaEventRecord(t_start);
    cudaMemcpy(d_params, h_params, 7*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_triangle_list, h_triangle_list, tl_size*sizeof(Triangle), cudaMemcpyHostToDevice);
    cudaMemcpy(d_light_positions, light_positions, n_lights*sizeof(Vec3f), cudaMemcpyHostToDevice);

    cudaEventRecord(k_start);
    trace_kernel <<< dimGrid, dimBlock >>> (d_params, d_triangle_list, d_light_positions, d_image);
    cudaEventRecord(k_stop);

    cudaMemcpy(image, d_image, width*height*sizeof(Vec3f), cudaMemcpyDeviceToHost);
    cudaEventRecord(t_stop);

    //Print timing data
    cudaEventSynchronize(k_stop);
    float k_time = 0;
    cudaEventElapsedTime(&k_time, k_start, k_stop);

    cudaEventSynchronize(t_stop);
    float t_time = 0;
    cudaEventElapsedTime(&t_time, t_start, t_stop);

    printf("GPU kernel time(ms): %f\n", k_time);
    printf("GPU total time(ms): %f\n", t_time);

    cudaFree(d_image);
    cudaFree(d_triangle_list);
    cudaFree(d_params);

    //Parallel program ends*/

    /*// Serial program begins
    cudaEventRecord(c_start);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -2);
            raydir.normalize();
            *pixel = trace(Vec3f(0,0,-7), raydir, h_triangle_list, tl_size);
        }
        std::cout << y << "\n";
    }
    cudaEventRecord(c_stop);

    //Print timing data
    cudaEventSynchronize(c_stop);
    float c_time = 0;
    cudaEventElapsedTime(&c_time, c_start, c_stop);

    printf("CPU total time(ms): %f\n", c_time);

    //Serial program ends*/

    //Write output to ppm file
    std::ofstream ofs(f_output, std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x/255)*255 ) <<
               (unsigned char)(std::min(float(1), image[i].y/255)*255 ) <<
               (unsigned char)(std::min(float(1), image[i].z/255)*255 );
    }
    ofs.close();

    //Free memory
    delete [] image;
    free(h_triangle_list);

}
void load_mesh_into_world(char const *, std::vector<Triangle *> &, Vec3f offset=Vec3f(0, 0, 0), const bool reflective=false);

int main(int argc, char const *argv[]){

    if (argc != 3){
        printf("Usage: ./raytacer <modelfile>.obj <output>.ppm\n");
        exit(0);
    }

    std::vector<Triangle *> triangle_list;
    load_mesh_into_world(argv[1], triangle_list, Vec3f(0, 0, 0));
    // load_mesh_into_world(argv[1], triangle_list, Vec3f(-1.5, 0, 0));
    // load_mesh_into_world(argv[1], triangle_list, Vec3f(1.5, 0, 0));


    // adding a ground plane
    load_mesh_into_world("plane.obj", triangle_list, Vec3f(0, 0, 0));

    // point light position coordinates
    Vec3f light_positions[] = {
        Vec3f(7, 7, -2),
        Vec3f(-3, -3, -5),
        Vec3f(10, 0, 0),
	Vec3f(-10, 0, 7)
    };

    // only works because light_positions is on stack and known at compile time
    int n_lights = sizeof(light_positions)/sizeof(Vec3f);

    printf("Rendering %d triangles...\n", triangle_list.size());
    render(triangle_list, light_positions, n_lights, argv[2]);

    return 0;
}


double det(double a1, double a2, double a3,
            double b1, double b2, double b3,
            double c1, double c2, double c3)
{
    double t1 = a1 * (b2*c3 - b3*c2);
    double t2 = a2 * (b1*c3 - b3*c1);
    double t3 = a3 * (b1*c2 - b2*c1);
    return t1 - t2 + t3;
}

int clamp(int what, int low, int high)
{
    if (what < low) return low;
    if (what > high) return high;
    return what;
}

void load_mesh_into_world(char const *mesh_filename, std::vector<Triangle *> &triangle_list, Vec3f offset, const bool reflective){

    std::ifstream objinfile(mesh_filename);

    std::string line;
    std::vector<Vec3f*> vertices;
    std::vector<Vec3f*> vertex_textures;
    Triangle* triangle;

    //Reading obj file
    while(getline(objinfile, line)){

        std::istringstream iss(line);
        std::string type_;
        iss >> type_;
        std::string fv1, fv2, fv3;

    //Create list of vertices
        if (type_.compare("v") == 0){

            double a, b, c;
            iss >> a >> b >> c;
            vertices.push_back(new Vec3f(a + offset.x, b + offset.y, c + offset.z));
        }
    //Create list of vertex textures
        else if (type_.compare("vt") == 0){
            double a, b;
            iss >> a >> b;
            vertex_textures.push_back(new Vec3f(a, b, 0));
        }
    //Create list of triangles
        else if (type_.compare("f") == 0){

            iss >> fv1 >> fv2 >> fv3;
            std::stringstream ssfv1(fv1);
            std::stringstream ssfv2(fv2);
            std::stringstream ssfv3(fv3);

            int v1, v2, v3;
            int vt1, vt2, vt3;
            ssfv1 >> v1;
            ssfv1.ignore();
            ssfv1 >> vt1;

            ssfv2 >> v2;
            ssfv2.ignore();
            ssfv2 >> vt2;

            ssfv3 >> v3;
            ssfv3.ignore();
            ssfv3 >> vt3;

            Vec3f *vertex1 = vertices[v1 - 1];
            Vec3f *vertex2 = vertices[v2 - 1];
            Vec3f *vertex3 = vertices[v3 - 1];

            triangle = new Triangle(*vertices[v1-1], *vertices[v2-1], *vertices[v3-1],
                                    *vertex_textures[vt1-1], *vertex_textures[vt2-1], *vertex_textures[vt3-1],
                                    reflective);

            triangle_list.push_back(triangle);
        }
    }
}
