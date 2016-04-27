#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>

#define MYASSERT(condition) if (!(condition)) { return; }


#define BLOCK_SIZE 32
#define HD __host__ __device__
#define WIDTH 512

HD float max_(float a, float b){ return (a < b) ? b : a; }
HD float min_(float a, float b){ return (a > b) ? b : a ;}
HD int clamp(int what, int low, int high);

#include "geometry.cuh"
#include "grid.cuh"

HD Vec3f trace(Ray ray, Triangle* triangle_list, int tl_size);
HD Vec3f fast_trace(Ray& ray, GridAccel* newGridAccel, int isDebugThread);

__global__
void trace_kernel (float* params, Vec3f* image,GridAccel* d_newGridAccel, Triangle* triangle_list, Vec3f _u, Vec3f _v, Vec3f _w, Vec3f camerapos){

    //Pixel coordinates
    int x = threadIdx.x + blockIdx.x*blockDim.x;
    int y = threadIdx.y + blockIdx.y*blockDim.y;

    //Unpack parameters
    int tl_size = (int) params[0];
	float width = params[1], height = params[2];
	float focal_distance = params[3];
	float aspectratio = params[4];
	Vec3f u = _u;
	Vec3f v = _v;
	Vec3f w = _w;

    Vec3f dir(0);
    dir = dir.add(w.negate().scale(focal_distance));
    float xw = aspectratio*(x - width/2.0 + 0.5)/width;
    float yw = (y - height/2.0 + 0.5)/height;
    dir = dir.add(u.scale(xw));
    dir = dir.add(v.scale(yw));
    dir.normalize();

    Ray ray(camerapos, dir, 0);

    //http://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-ray-tracing/ray-tracing-practical-example

    // trace the ray
    // image[y*WIDTH + x] = trace(ray, triangle_list, tl_size);
    image[y*WIDTH + x] = fast_trace(ray, d_newGridAccel, x == 275  && y == 240);

}
 
Vec3f trace(Ray ray, Triangle* triangle_list, int tl_size)
{
    // TODO: pass as pointer in parameters
    // base color, kd, ks, spec_alpha, ka, reflective
    material materials[2] = {};
    init_material(&materials[0], Vec3f(255, 0, 0), 10.0f, 10.0f, 1.25, 0.3f, false); // blub
    init_material(&materials[1], Vec3f(0, 0, 255), 1.0f, 1.5f, 1.25, 0.3f, false); // plane

    Intersection* isect;
    Vec3f rayorig = ray.orig;
    Vec3f raydir = ray.raydir;

    double tnear = INFINITY;
    Vec3f normal;
 
    Triangle triangle_near(Vec3f(100), Vec3f(100), Vec3f(100),Vec3f(100),Vec3f(100), Vec3f(100), 0);
    bool hit = false;
    for (unsigned int i = 0; i < tl_size; ++i) {
        hit |= triangle_list[i].Intersect(ray, isect, triangle_near, tnear, normal);
    }
    if (!hit)
        return Vec3f(0);

    Vec3f light_pos(2, 5, 0);
    // light_pos = ray.orig;
    // light_pos.z = -1*light_pos.z;

    Vec3f poi = rayorig.add( raydir.scale(tnear) );
    Vec3f v = raydir.negate().normalize();
    Vec3f l = light_pos.subtract(poi).normalize();
    Vec3f h = v.add(l).normalize();

    // Simple blinn phong shading
 
    material mat = materials[triangle_near.material_index];
    Vec3f base_color = mat.base_color;
    float kd = mat.kd;
    float ks = mat.ks;
    float spec_alpha = mat.spec_alpha;
    float ka = mat.ka;

    Vec3f diffuse = base_color.multiply(max_(float(0),normal.dotProduct(l))).scale(kd);
    Vec3f specular = base_color.multiply(pow(max_(float(0), normal.dotProduct(h)),spec_alpha)).scale(ks);
    Vec3f ambient = base_color.scale(ka);
    Vec3f color = diffuse.add(specular).add(ambient);
    
    Vec3f shadow_ray_dir = light_pos.subtract(poi).normalize();

    Ray shadow_ray(poi, shadow_ray_dir, eps + 2);
    bool in_shadow = false;
    {
        Triangle temp_tri_near(Vec3f(100), Vec3f(100), Vec3f(100),Vec3f(100),Vec3f(100), Vec3f(100), 0);
        double temp_tnear = INFINITY;
        Vec3f temp_normal(0, 0, 0);
        for (unsigned int i = 0; i < tl_size; ++i) {
            if(triangle_list[i].Intersect(shadow_ray, isect, temp_tri_near, temp_tnear, temp_normal)){
                in_shadow = true;
                break;
            }
        }
    }
    if (in_shadow){
        return color.scale(0.5f);
    }
    else{

    }

    return color;
}


Vec3f fast_trace(Ray& ray, GridAccel* newGridAccel, int isDebugThread){

    // TODO: pass as pointer in parameters
    // base color, kd, ks, spec_alpha, ka, reflective
    material materials[2] = {};
    init_material(&materials[0], Vec3f(255, 0, 0), 10.0f, 10.0f, 1.25, 0.3f, false); // blub
    init_material(&materials[1], Vec3f(0, 0, 255), 1.0f, 1.5f, 1.25, 0.3f, false); // plane

    Intersection* isect = NULL;
    Vec3f rayorig = ray.orig, raydir = ray.raydir;

    //Practically infinity
    Triangle triangle_near(Vec3f(100), Vec3f(100), Vec3f(100),Vec3f(100),Vec3f(100), Vec3f(0), 0);
    double tnear = INFINITY;
    Vec3f normal(0);

    bool hitSomething = newGridAccel->Intersect(ray, isect, triangle_near, tnear, normal, isDebugThread);

    if (!hitSomething)
	    return Vec3f(0);


    // assume only 1 light over here.
    Vec3f light_pos(2, 5, 0);
    // light_pos = ray.orig;
    // light_pos.z = -1*light_pos.z;
 
    Vec3f poi = rayorig.add( raydir.scale(tnear) );
    Vec3f v = raydir.negate().normalize();
    Vec3f l = light_pos.subtract(poi).normalize();
    Vec3f h = v.add(l).normalize();
 
    // Simple blinn phong shading
    material mat = materials[triangle_near.material_index];
    Vec3f base_color = mat.base_color;
    float kd = mat.kd;
    float ks = mat.ks;
    float spec_alpha = mat.spec_alpha;
    float ka = mat.ka;

    Vec3f diffuse = base_color.multiply(max_(float(0),normal.dotProduct(l))).scale(kd);
    Vec3f specular = base_color.multiply(pow(max_(float(0), normal.dotProduct(h)),spec_alpha)).scale(ks);
    Vec3f ambient = base_color.scale(ka);
    Vec3f color = diffuse.add(specular).add(ambient);

    Vec3f shadow_ray_dir = light_pos.subtract(poi).normalize();

    Ray shadow_ray(poi, shadow_ray_dir, eps + 0.02f);
    
    bool in_shadow = false;
    {
        Triangle temp_tri_near(Vec3f(100), Vec3f(100), Vec3f(100),Vec3f(100),Vec3f(100), Vec3f(0), 0);
        double temp_tnear = INFINITY;
        Vec3f temp_normal(0);

        in_shadow = newGridAccel->Intersect(shadow_ray, isect, temp_tri_near, temp_tnear, temp_normal, isDebugThread);
    }
    if (in_shadow){
        return color.scale(0.5f);
    }
    

    return color;
}

void render(std::vector<Triangle*> &triangle_list){

    //Define image size, calculate camera view parameters
    int width = WIDTH, height = WIDTH;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    Vec3f camera_pos(4+4+10, 4+4+10, 5+4+10);  //camera_pos(0, -500, -100); nefertiti
    Vec3f camera_target(0, 0, 0);
    Vec3f camera_up(0, -1, 0);
    float fov = 60;

    camera_up.normalize();
    Vec3f line_of_sight = camera_target.subtract(camera_pos);
    Vec3f w = line_of_sight.negate().normalize();
    Vec3f u = camera_up.crossProduct(w).normalize();
    Vec3f v = w.crossProduct(u).normalize();
    float focal_height = 1.0f;
    float aspectratio = float(width)/float(height);
    float focal_width = focal_height * aspectratio;
    float focal_distance = focal_height/(2.0 * tan(fov * M_PI/(180.0 * 2.0)));

    int tl_size = triangle_list.size();

    //Cuda events for timing data
    cudaEvent_t c_start, c_stop, t_start, t_stop, k_start, k_stop;
    cudaEventCreate(&c_start);
    cudaEventCreate(&c_stop);
    cudaEventCreate(&t_start);
    cudaEventCreate(&t_stop);
    cudaEventCreate(&k_start);
    cudaEventCreate(&k_stop);

	int params_size = 5; //Used to control memory allocation and memcpy down below. Should match size of h_params
    float h_params[5] = {tl_size, width, height, focal_distance, aspectratio};

    Triangle* h_triangle_list = (Triangle*)malloc(tl_size*sizeof(Triangle));
    for(int i = 0;i<tl_size;i++)
        h_triangle_list[i] = *triangle_list[i];

    GridAccel* newGridAccel = new GridAccel(h_triangle_list, tl_size);

	int totalVoxels = 1;
	for(int i=0;i<3;i++)
		totalVoxels *= newGridAccel->nVoxels[i];

    Voxel** h_voxels = (Voxel**)malloc(sizeof(Voxel*)*totalVoxels);

    //Parallel program begins
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(width/dimBlock.x, height/dimBlock.y);

    float* d_params;
    Triangle* d_triangle_list;
    Vec3f *d_image;
	GridAccel* d_newGridAccel;
	Voxel** d_voxels;

    //Copy parameters needed to set up camera view, triangle list
    cudaMalloc(&d_params, params_size*sizeof(float));
    cudaMalloc(&d_triangle_list, tl_size*sizeof(Triangle));
    cudaMalloc(&d_image, width*height*sizeof(Vec3f));
	cudaMalloc(&d_voxels, totalVoxels*sizeof(Voxel*));

	// Copying voxels to device memory

    cudaEventRecord(t_start);
	int cnt = 0;

	for(int i = 0;i<totalVoxels;i++){

        if(newGridAccel->voxels[i] != NULL){

			Triangle* d_voxel_triangle_list;
           	cudaMalloc(&d_voxel_triangle_list, newGridAccel->voxels[i]->voxelListSize*sizeof(Triangle));
			cudaMemcpy(d_voxel_triangle_list, newGridAccel->voxels[i]->triangleList, newGridAccel->voxels[i]->voxelListSize*sizeof(Triangle), cudaMemcpyHostToDevice);

			h_voxels[i] = (Voxel*)malloc(sizeof(Voxel));
			h_voxels[i]->voxelListSize = newGridAccel->voxels[i]->voxelListSize;
			h_voxels[i]->triangleList = d_voxel_triangle_list;

			Voxel* d_voxel_elem;
			cudaMalloc(&d_voxel_elem, sizeof(Voxel));
			cudaMemcpy(d_voxel_elem, h_voxels[i], sizeof(Voxel), cudaMemcpyHostToDevice);

			h_voxels[i] = d_voxel_elem;

        }
		cnt++;
    }

    cudaMemcpy(d_voxels, h_voxels, totalVoxels*sizeof(Voxel*), cudaMemcpyHostToDevice);

    //Allocate memory for GridAccel
    cudaMalloc((void **)&d_newGridAccel, sizeof(GridAccel));
    newGridAccel->triangleList = d_triangle_list;
    newGridAccel->voxels = d_voxels;

    cudaMemcpy(d_params, h_params, params_size*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_triangle_list, h_triangle_list, tl_size*sizeof(Triangle), cudaMemcpyHostToDevice);

	//Copy GridAccel
	cudaMemcpy(d_newGridAccel, newGridAccel, sizeof(GridAccel), cudaMemcpyHostToDevice);

	cudaEventRecord(k_start);
    trace_kernel <<< dimGrid, dimBlock >>> (d_params, d_image, d_newGridAccel, d_triangle_list, u, v , w, camera_pos);
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
	cudaFree(d_voxels);
	cudaFree(d_newGridAccel);

    //Parallel program ends*/

    /*// Serial program begins : NB : Camera has changed. This will not work now.
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
    std::ofstream ofs("./nefertiti0.ppm", std::ios::out | std::ios::binary);
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
	free(h_voxels);

}


int main(){

    std::vector<Triangle*> triangle_list;

    // params: filename, &triangle_list, format_has_vt, material_index, offset, scale

    // load_mesh("nefertiti_triangulated.obj", triangle_list, false, Vec3f(255, 0, 0), false, Vec3f(0, 0, 0)); // Incorrect call.
    // load_mesh("blub_triangulated.obj", triangle_list, true, 0, Vec3f(0, 0, 0), 5);
    load_mesh("blub_triangulated.obj", triangle_list, true, 0, Vec3f(0, 0, 0), 5);
    load_mesh("spot_triangulated.obj", triangle_list, true, 0, Vec3f(-2, 0, 0), 5);
    load_mesh("blub_triangulated.obj", triangle_list, true, 0, Vec3f(2, 0, 0), 5);
    load_mesh("plane.obj", triangle_list, true, 1, Vec3f(0, 0.4, 0), 3);

    std::cout << "Rendering " << triangle_list.size() << " triangles" << std::endl;
    render(triangle_list);

    return 0;
}

HD double det(double a1, double a2, double a3,
            double b1, double b2, double b3,
            double c1, double c2, double c3)
{
    double t1 = a1 * (b2*c3 - b3*c2);
    double t2 = a2 * (b1*c3 - b3*c1);
    double t3 = a3 * (b1*c2 - b2*c1);
    return t1 - t2 + t3;
}

HD int clamp(int what, int low, int high)
{
    if (what < low) return low;
    if (what > high) return high;
    return what;
}

void load_mesh(const char * filename, std::vector<Triangle *> &triangle_list, bool format_has_vt, int material_index, Vec3f offset, int scale)
{
    std::ifstream objinfile(filename);

    std::string line;
    std::vector<Vec3f*> vertices;
    std::vector<Vec3f*> vertex_textures;
    Triangle* triangle;

    while(getline(objinfile, line)){

        std::istringstream iss(line);
        std::string type_;
        iss >> type_;
        std::string fv1, fv2, fv3;
        if (type_.compare("v") == 0){

            double a, b, c;
            iss >> a >> b >> c;
            Vec3f *v_new = new Vec3f(scale*(a + offset.x), scale*(b + offset.y), scale*(c + offset.z));
            vertices.push_back(v_new);
        }
        else if (type_.compare("vt") == 0){
            double a, b;
            iss >> a >> b;
            vertex_textures.push_back(new Vec3f(a, b, 0));
        }
        else if (type_.compare("f") == 0){

            iss >> fv1 >> fv2 >> fv3;
            std::stringstream ssfv1(fv1);
            std::stringstream ssfv2(fv2);
            std::stringstream ssfv3(fv3);

            int v1, v2, v3;
            int vt1, vt2, vt3;
            ssfv1 >> v1;

            ssfv2 >> v2;

            ssfv3 >> v3;

            Vec3f *vertex1 = vertices[v1 - 1];
            Vec3f *vertex2 = vertices[v2 - 1];
            Vec3f *vertex3 = vertices[v3 - 1];

            if (format_has_vt){
                ssfv1.ignore();
                ssfv1 >> vt1;
                ssfv2.ignore();
                ssfv2 >> vt2;
                ssfv3.ignore();
                ssfv3 >> vt3;

                triangle = new Triangle(*vertices[v1-1], *vertices[v2-1], *vertices[v3-1],
                                        *vertex_textures[vt1-1], *vertex_textures[vt2-1], *vertex_textures[vt3-1],
                                        material_index);
            } else{
                triangle = new Triangle(*vertices[v1-1], *vertices[v2-1], *vertices[v3-1],
                                        0, 0, 0,
                                        material_index);
            }

            triangle_list.push_back(triangle);
        }
    }

    objinfile.close();
}
