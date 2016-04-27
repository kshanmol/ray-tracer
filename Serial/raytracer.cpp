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

int clamp(int what, int low, int high);

#include "geometry.h"
#include "grid.h"

material materials[] = {
    // base color, kd, ks, spec_alpha, ka, reflective
    material{ Vec3f(255, 0, 0), 10.0f, 10.0f, 1.25, 0.3f, false}, // blub
    material{ Vec3f(0, 0, 255), 1.0f, 1.5f, 1.25, 0.3f, false} // plane
};

Vec3f reflect(const Vec3f &I, const Vec3f &N){
        return I.subtract(N.scale(2*I.dotProduct(N))).negate();
}

Vec3f trace(Ray ray, const std::vector<Triangle*> &triangle_list)
{
    Intersection* isect = new Intersection();
    Vec3f rayorig = ray.orig;
    Vec3f raydir = ray.raydir;

    double tnear = INFINITY;
    Vec3f normal;

    Triangle triangle_near(Vec3f(100), Vec3f(100), Vec3f(100),Vec3f(100),Vec3f(100), Vec3f(100), 0);
    bool hit = false;
    for (unsigned int i = 0; i < triangle_list.size(); ++i) {
        hit |= triangle_list[i]->Intersect(ray, isect, triangle_near, tnear, normal);
    }
    if (!hit)
        return Vec3f(0);

    /*
    Vector3D v = -incident.getDirection();
    v.normalize();
    Vector3D l = world->getLightSourcePosition() - incident.getPosition();
    l.normalize();
    Vector3D h = (v+l);
    h.normalize();
    Vector3D normal = incident.getNormal();
    normal.normalize();

    Color I = world->getLightSourceIntensity();
    Color L = I*fmax(0,dotProduct(normal,l))*kd + I*pow(fmax(0, dotProduct(normal,h)),this->n)*ks + I*ka;
    */

    Vec3f light_pos(5, -5, 0);
    light_pos = ray.orig;
    // std::cout << ray.orig << std::endl;

    Vec3f poi = rayorig.add( raydir.scale(tnear) );
    Vec3f v = raydir.negate().normalize();
    Vec3f l = light_pos.subtract(poi).normalize();
    Vec3f h = v.add(l).normalize();

    normal = normal.negate();

    // Simple blinn phong shading

    material mat = materials[triangle_near.material_index];
    Vec3f base_color = mat.base_color;
    float kd = mat.kd;
    float ks = mat.ks;
    float spec_alpha = mat.spec_alpha;
    float ka = mat.ka;

    Vec3f diffuse = base_color.multiply(std::max(float(0),normal.dotProduct(l))).scale(kd);
    Vec3f specular = base_color.multiply(pow(std::max(float(0), normal.dotProduct(h)),spec_alpha)).scale(ks);
    // std::cout << normal.dotProduct(h) << std::endl;
    Vec3f ambient = base_color.scale(ka);
    return diffuse.add(specular).add(ambient);

    // Vec3f poi = rayorig.add( raydir.scale(tnear) );
    // Vec3f eye = rayorig.subtract(poi).normalize();  //raydir.negate();
    // Vec3f l = light_pos.subtract(poi).normalize();
    // Vec3f half = eye.add(l).normalize();
    // Vec3f n = triangle_near->getNormal(poi).normalize();


    // Vec3f diffuse = color.scale(kd * std::max(float(0), n.dotProduct(l.normalize())));
    // Vec3f specular = color.scale(ks * pow(std::max(float(0), reflect(l,n).dotProduct(raydir.negate())), spec_alpha));
    // Vec3f ambient = Vec3f(40.0f);


    // // actual
    // return diffuse.add(specular).add(ambient);

}

Vec3f fast_trace(Ray& ray, GridAccel* newGridAccel)
{
    /*
	Intersection* isect = new Intersection();
	Vec3f rayorig = ray.orig, raydir = ray.raydir;
	global_t = INFINITY;

    Vec3f normal
	bool hitSomething = newGridAccel->Intersect(ray, isect);
    if (!hitSomething)
        return Vec3f(0); // background colour.

    // Simple blinn phong shading
    Vec3f color = global_triangle_near->material;
    float kd = 2.0f;
    float ks = 50.0f;
    float ka = 0.2f;
    float spec_alpha = 4;
    Vec3f light_pos(5, -5, 2);
    float light_intensity = 255;
    float shadow_scale = 0.2f;

    // std::cout << global_t << std::endl;

    // blinn phong
    Vec3f poi = rayorig.add( raydir.scale(global_t) );
    Vec3f v = raydir.negate().normalize();
    Vec3f l = light_pos.subtract(poi).normalize();
    Vec3f h = v.add(l).normalize();
    Vec3f normal = global_triangle_near->getNormalMod();

    Vec3f diffuse = color.scale(kd * std::max(float(0), normal.dotProduct(l))).scale(light_intensity);
    Vec3f specular = color.scale(ks * pow(std::max(float(0), normal.dotProduct(h)), spec_alpha)).scale(light_intensity);
    Vec3f ambient = color.scale(ka);

    color = specular.add(diffuse).add(ambient);

    // shadow
    Vec3f shadow_dir = light_pos.subtract(poi).normalize();
    Ray shadow_ray(poi, shadow_dir, eps);

    isect->use_eps = true;
    hitSomething = newGridAccel->Intersect(shadow_ray, isect);
    // if (hitSomething)
    //     return Vec3f(0, 255, 0);
        // color = color.scale(shadow_scale);
    // isect->use_eps = false;


    // assume only 1 light over here.
    return color;
    */
    return Vec3f(255, 0, 0);
}

void render(const std::vector<Triangle*> &triangle_list){

    GridAccel* newGridAccel = new GridAccel(triangle_list);

    Vec3f camera_pos(4+2+2, 4+2+2, 5+2+2);
    Vec3f camera_target(0, 0, 0);
    Vec3f camera_up(0, -1, 0);
    float fov = 60;
    int width = 256, height = 256;

    camera_up.normalize();
    Vec3f line_of_sight = camera_target.subtract(camera_pos);
    Vec3f w = line_of_sight.negate().normalize();
    Vec3f u = camera_up.crossProduct(w).normalize();
    Vec3f v = w.crossProduct(u).normalize();
    float focal_height = 1.0f;
    float aspectratio = float(width)/float(height);
    float focal_width = focal_height * aspectratio;
    float focal_distance = focal_height/(2.0 * tan(fov * M_PI/(180.0 * 2.0)));


    Vec3f *image = new Vec3f[width * height], *pixel = image;

    // Trace rays
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x, ++pixel) {

            Vec3f dir(0);
            dir = dir.add(w.negate().scale(focal_distance));
            float xw = aspectratio*(x - width/2.0 + 0.5)/width;
            float yw = (y - height/2.0 + 0.5)/height;
            dir = dir.add(u.scale(xw));
            dir = dir.add(v.scale(yw));
            dir.normalize();

            Ray ray(camera_pos, dir, 0);

            // *pixel = trace(ray, Vec3f(0,0.1,-7), raydir, triangle_list);
            // *pixel = fast_trace(ray, newGridAccel);
            *pixel = trace(ray, triangle_list);
        }
        std::cout << y << "\n";
    }

    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./test.ppm", std::ios::out | std::ios::binary);
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

    std::vector<Triangle*> triangle_list;

    //params: filename, &triangle_list, format_has_vt, material_index, offset, scale
    load_mesh("blub_triangulated.obj", triangle_list, true, 0, Vec3f(0, 0, 0), 5);
    load_mesh("plane.obj.new", triangle_list, true, 1, Vec3f(0, 0.85, 0), 10);

    std::cout << "Rendering " << triangle_list.size() << " triangles" << std::endl;
    render(triangle_list);

    return 0;
}


double det(double data_3x3[3][3])
{
    return data_3x3[0][0]*((data_3x3[1][1]*data_3x3[2][2]) - (data_3x3[2][1]*data_3x3[1][2])) -data_3x3[0][1]*(data_3x3[1][0]*data_3x3[2][2] - data_3x3[2][0]*data_3x3[1][2]) + data_3x3[0][2]*(data_3x3[1][0]*data_3x3[2][1] - data_3x3[2][0]*data_3x3[1][1]);
}

int clamp(int what, int low, int high)
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
            Vec3f *v_new = new Vec3f(scale*(a + offset.x), scale*(b + offset.y), scale*(c + offset.z)); // TODO: check.
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
