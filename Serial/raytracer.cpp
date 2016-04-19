#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>

static const float eps = 1e-8;
int clamp(int what, int low, int high);

#include "geometry.h"
#include "grid.h"

struct Intersection{


};

typedef struct Intersection Intersection;

Vec3f reflect(const Vec3f &I, const Vec3f &N){
        return I.subtract(N.scale(2*I.dotProduct(N))).negate();
}

Vec3f trace(Ray ray, Vec3f rayorig, Vec3f raydir,
            const std::vector<Triangle*> &triangle_list)
{
    float tnear = INFINITY,
          beta = INFINITY,
          gamma = INFINITY;

    const Triangle* triangle_near = NULL;
    for (unsigned int i = 0; i < triangle_list.size(); ++i) {
        float t0 = INFINITY,
              beta_ = INFINITY,
              gamma_ = INFINITY;
        if (triangle_list[i]->rayTriangleIntersect(rayorig, raydir, t0, beta_, gamma_)) {
            if (t0 < tnear) {
                tnear = t0;
                triangle_near = triangle_list[i];
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
    Vec3f light_pos(5, 5, -2);

    Vec3f poi = rayorig.add( raydir.scale(tnear) );
    Vec3f eye = rayorig.subtract(poi).normalize();  //raydir.negate();
    Vec3f l = light_pos.subtract(poi).normalize();
    Vec3f half = eye.add(l).normalize();
    Vec3f n = triangle_near->getNormal(poi).normalize();

   
    Vec3f diffuse = color.scale(kd * std::max(float(0), n.dotProduct(l.normalize())));
    Vec3f specular = color.scale(ks * pow(std::max(float(0), reflect(l,n).dotProduct(raydir.negate())), spec_alpha));
    Vec3f ambient = Vec3f(40.0f);


    // actual
    return diffuse.add(specular).add(ambient);

}

Vec3f fast_trace(Ray& ray, GridAccel* newGridAccel, int isDebugThread)
{

	Intersection* isect;
	Vec3f rayorig = ray.orig, raydir = ray.raydir;
	global_t = INFINITY;

	bool hitSomething = newGridAccel->Intersect(ray, isect, isDebugThread);    
    if (!hitSomething)
        return Vec3f(0);

    // Simple blinn phong shading
    Vec3f color(200.0);
    float kd = 0.3f;
    float ks = 0.5f;
    float spec_alpha = 4;

    // assume only 1 light over here.
    Vec3f light_pos(5, 5, -2);

    Vec3f poi = rayorig.add( raydir.scale(global_t) );
    Vec3f eye = rayorig.subtract(poi).normalize();  //raydir.negate();
    Vec3f l = light_pos.subtract(poi).normalize();
    Vec3f half = eye.add(l).normalize();
    Vec3f n = global_triangle_near->getNormal(poi).normalize();

   
    Vec3f diffuse = color.scale(kd * std::max(float(0), n.dotProduct(l.normalize())));
    Vec3f specular = color.scale(ks * pow(std::max(float(0), reflect(l,n).dotProduct(raydir.negate())), spec_alpha));
    Vec3f ambient = Vec3f(40.0f);


    // actual
    return diffuse.add(specular).add(ambient);

}

void render(const std::vector<Triangle*> &triangle_list){

    GridAccel* newGridAccel = new GridAccel(triangle_list);

    int width = 64, height = 64;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 18);

    // Trace rays
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, 2);
            raydir.normalize();
            Ray ray(Vec3f(0,0,-7), raydir, 0);
            //*pixel = trace(ray, Vec3f(0,0.1,-7), raydir, triangle_list);
            *pixel = fast_trace(ray, newGridAccel, x == 32 && y == 32);
        }
        //std::cout << y << "\n";
    }
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./test1.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x/255)*255 ) <<
               (unsigned char)(std::min(float(1), image[i].y/255)*255 ) <<
               (unsigned char)(std::min(float(1), image[i].z/255)*255 );
    }
    ofs.close();
    delete [] image;
}


int main(){

    std::ifstream objinfile("spot_triangulated.obj");

    std::string line;
    std::vector<Vec3f*> vertices;
    std::vector<Vec3f*> vertex_textures;
    std::vector<Triangle*> triangle_list;
    Triangle* triangle;

    while(getline(objinfile, line)){

        std::istringstream iss(line);
        std::string type_;
        iss >> type_;
        std::string fv1, fv2, fv3;
        if (type_.compare("v") == 0){

            double a, b, c;
            iss >> a >> b >> c;
            vertices.push_back(new Vec3f(a, b, c));
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
