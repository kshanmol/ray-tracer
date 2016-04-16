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

static const float eps = 1e-8;
int clamp(int what, int low, int high);

#include "geometry.cuh"
#include "grid.cuh"

HD Vec3f trace(Vec3f rayorig, Vec3f raydir, Triangle* triangle_list, int tl_size, GridAccel* d_newGridAccel);

HD float max_(float a, float b){ return (a < b) ? b : a; }

__global__
void trace_kernel (float* params, Triangle* triangle_list, Vec3f* image,GridAccel* d_newGridAccel){

    //Pixel coordinates
    int x = threadIdx.x + blockIdx.x*blockDim.x;
    int y = threadIdx.y + blockIdx.y*blockDim.y;	

    //Unpack parameters
    float invWidth = params[0], invHeight = params[1];
    float fov = params[2], aspectratio = params[3];
    float angle = params[4];
    int tl_size = (int) params[5];

    //http://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-ray-tracing/ray-tracing-practical-example

    //Set up camera view
    float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
    float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
    Vec3f raydir(xx, yy, -2);
    raydir.normalize();

    //Trace ray
    image[y*WIDTH + x] = trace(Vec3f(0,0,-7), raydir, triangle_list, tl_size, d_newGridAccel);

}

Vec3f trace(Vec3f rayorig, Vec3f raydir, Triangle* triangle_list, int tl_size, GridAccel* newGridAccel)
{

    //Ray triangle intersection
    float tnear = INFINITY,
          beta = INFINITY,
          gamma = INFINITY;

    const Triangle* triangle_near = NULL;
    for (unsigned int i = 0; i < tl_size; ++i) {
        float t0 = INFINITY,
              beta_ = INFINITY,
              gamma_ = INFINITY;
        if (triangle_list[i].rayTriangleIntersect(rayorig, raydir, t0, beta_, gamma_)) {
            if (t0 < tnear) {
                tnear = t0;
                triangle_near = &triangle_list[i];
                beta = beta_;
                gamma = gamma_;
            }
        }
    }
    if (!triangle_near)
        return Vec3f(0);
 
    // Simple blinn phong shading
    Vec3f color(200.0);
    float kd = 0.3f;
    float ks = 0.5f;
    float spec_alpha = 4;

    // assume only 1 light over here.
    Vec3f light_pos(7, 7, -2);

    Vec3f poi = rayorig.add( raydir.scale(tnear) );
    Vec3f eye = rayorig.subtract(poi).normalize();  //raydir.negate();
    Vec3f l = poi.subtract(light_pos).normalize();
    Vec3f half = eye.add(l).normalize();
    Vec3f n = triangle_near->getNormal(poi).normalize();

    Vec3f diffuse = color.scale(kd * max_(float(0), n.dotProduct(l.normalize())));
    Vec3f specular = color.scale(ks * pow(max_(float(0), reflect(l,n).dotProduct(raydir.negate())), spec_alpha));
    Vec3f ambient = Vec3f(40.0);

    //return specular;
    // actual
    return diffuse.add(specular).add(ambient);

}

void render(std::vector<Triangle*> &triangle_list){

	//GridAccel* newGridAccel = new GridAccel(triangle_list);

    //Define image size, calculate camera view parameters
    int width = 1024, height = 1024;
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
   
    float h_params[6] = {invWidth, invHeight, fov, aspectratio, angle, tl_size};

    Triangle* h_triangle_list = (Triangle*)malloc(tl_size*sizeof(Triangle));
    for(int i = 0;i<tl_size;i++)
        h_triangle_list[i] = *triangle_list[i];

    GridAccel* newGridAccel = new GridAccel(h_triangle_list, tl_size);
	GridAccel* h_newGridAccel = new GridAccel(h_triangle_list, tl_size);	

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
    cudaMalloc(&d_params, 6*sizeof(float));    
    cudaMalloc(&d_triangle_list, tl_size*sizeof(Triangle));
    cudaMalloc(&d_image, width*height*sizeof(Vec3f));
	cudaMalloc(&d_voxels, totalVoxels*sizeof(Voxel*));

	
	// Uncommenting the following causes a seg fault on iteration number 878
	/*int cnt = 0;

	for(int i = 0;i<totalVoxels;i++){

        if(newGridAccel->voxels[i] != NULL){

			Voxel* d_voxel_elem;
            cudaMalloc(&d_voxel_elem, sizeof(Voxel));
			cudaMemcpy(d_voxel_elem, newGridAccel->voxels[i], sizeof(Voxel), cudaMemcpyHostToDevice);
		
			Triangle* voxel_triangle_list;
           	cudaMalloc(&voxel_triangle_list, newGridAccel->voxels[i]->voxelListSize*sizeof(Triangle));
			cudaMemcpy(voxel_triangle_list, newGridAccel->voxels[i]->triangleList, newGridAccel->voxels[i]->voxelListSize*sizeof(Triangle), cudaMemcpyHostToDevice); 
			
			h_voxels[i] = d_voxel_elem;
			h_voxels[i]->triangleList = (Triangle*)malloc(newGridAccel->voxels[i]->voxelListSize*sizeof(Triangle));
			h_voxels[i]->triangleList = voxel_triangle_list;

			cudaFree(d_voxel_elem);
			cudaFree(voxel_triangle_list);
        }
		cnt++;
    }*/

    cudaMemcpy(d_voxels, h_voxels, totalVoxels*sizeof(Voxel*), cudaMemcpyHostToDevice);

    //Allocate memory for GridAccel
    cudaMalloc((void **)&d_newGridAccel, sizeof(GridAccel));
    newGridAccel->triangleList = d_triangle_list;
    newGridAccel->voxels = d_voxels;

    cudaEventRecord(t_start);
    cudaMemcpy(d_params, h_params, 6*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_triangle_list, h_triangle_list, tl_size*sizeof(Triangle), cudaMemcpyHostToDevice);
	
	//Copy GridAccel
	cudaMemcpy(d_newGridAccel, newGridAccel, sizeof(GridAccel), cudaMemcpyHostToDevice);
   

	cudaEventRecord(k_start);
    trace_kernel <<< dimGrid, dimBlock >>> (d_params, d_triangle_list, d_image, d_newGridAccel);
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
    std::ofstream ofs("./gpu_trial2.ppm", std::ios::out | std::ios::binary);
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


int main(){

    std::ifstream objinfile("spot_triangulated.obj");

    std::string line;
    std::vector<Vec3f*> vertices;
    std::vector<Vec3f*> vertex_textures;
    std::vector<Triangle*> triangle_list;
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
            vertices.push_back(new Vec3f(a, b, c));
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
                                    *vertex_textures[vt1-1], *vertex_textures[vt2-1], *vertex_textures[vt3-1]);

            triangle_list.push_back(triangle);
        }
    }

    render(triangle_list);

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

