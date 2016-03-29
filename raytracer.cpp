#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>
#include <png++/png.hpp>

static const float eps = 1e-8;
double det(double a1, double a2, double a3,
            double b1, double b2, double b3,
            double c1, double c2, double c3);
int clamp(int what, int low, int high);

template<typename T>
class Vec3{
public:

    T x, y, z;
    //Constructors
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T val) : x(val), y(val), z(val) {}
    Vec3(T xval, T yval, T zval) : x(xval), y(yval), z(zval) {}

    Vec3& normalize(){
        T nor2 = length2();
        if (nor2 > 0) {
            T nor_inv = 1 / sqrt(nor2);
            x *= nor_inv, y *= nor_inv, z *= nor_inv;
        }
        return *this;
    }

    T dotProduct(const Vec3<T> &v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    Vec3<T> crossProduct(const Vec3<T> &v) const {

        T tmpX = y * v.z - z * v.y;
        T tmpY = z * v.x - x * v.z;
        T tmpZ = x * v.y - y * v.x;
        return Vec3<T>(tmpX, tmpY, tmpZ );
    }

    T length2(){
        return x * x + y * y + z * z;
    }
    T length(){
        return sqrt(length2());
    }

    Vec3<T> scale(const T &f) const {
        return Vec3<T>(x * f, y * f, z * f);
    }

    Vec3<T> multiply(const Vec3<T> &v) const {
        return Vec3<T>(x * v.x, y * v.y, z * v.z);
    }

    Vec3<T> subtract(const Vec3<T> &v) const {
        return Vec3<T>(x - v.x, y - v.y, z - v.z);
    }

    Vec3<T> add(const Vec3<T> &v) const {
        return Vec3<T>(x + v.x, y + v.y, z + v.z);
    }

    Vec3<T> negate() const {
        return Vec3<T>(-x, -y, -z);
    }

    //Helper function to format display
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v){
        os << "(" << v.x << " " << v.y << " " << v.z << ")";
        return os;
    }
};

typedef Vec3<float> Vec3f;

class Triangle{
public:

    Vec3f v0, v1, v2;
    Vec3f tv0, tv1, tv2; // texture coordinates of vertices

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

    //Not our stuff yet.
    bool rayTriangleIntersect(const Vec3f &orig, const Vec3f &dir, float &t, float &beta, float &gamma){

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

    Vec3f getNormal(Vec3f point) const
    {
        // from http://math.stackexchange.com/a/137551
        Vec3f p = point.subtract(v1);
        Vec3f q = v0.subtract(v2);

        // use point here. !!

        return Vec3f(
                    p.y * q.z - p.z * q.y,
                    -1 * (p.x*q.z - p.z * q.x),
                    p.x*q.y - p.y*q.x
                    );
    }
};

Vec3f trace(Vec3f rayorig, Vec3f raydir,
            const std::vector<Triangle*> &triangle_list,
            const png::image< png::rgb_pixel > texture_map)
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

    // map texture
    int u = (int) (((1 - beta - gamma)*triangle_near->tv1.x + beta*triangle_near->tv0.x + gamma*triangle_near->tv2.x) * texture_map.get_width());
    int v = (int) (((1 - beta - gamma)*triangle_near->tv1.y + beta*triangle_near->tv0.y + gamma*triangle_near->tv2.y) * texture_map.get_height());

    u = clamp(u, 0, texture_map.get_width() - 1);
    v = clamp(u, 0, texture_map.get_height() - 1);

    png::rgb_pixel texture = texture_map.get_pixel(u, v);

    // Simple blinn phong shading
    Vec3f color(fabs(100));
    float kd = 2;
    float ks = 4;
    float spec_alpha = 2;

    // assume only 1 light over here.
    Vec3f light_pos(1.5, 1.5, 1.5);

    Vec3f v = raydir.negate();
    Vec3f poi = rayorig.add( raydir.scale(tnear) );
    Vec3f l = light_pos.subtract(poi);
    Vec3f half = v.add(l).normalize();
    Vec3f n = triangle_near->getNormal(poi);

    // rayorig.negate();
    // raydir.negate();
    // Vec3f half = (rayorig.add(raydir)).normalize();
    // Vec3f normal = triangle_near->getNormal(rayorig.negate().add( raydir.negate().scale(tnear) ));

    Vec3f diffuse = color.scale(kd * std::max(float(0), n.dotProduct(l)));
    Vec3f specular = color.scale(ks * pow(std::max(float(0), n.dotProduct(h)), spec_alpha));

    // return color.scale(kd*std::max(float(0), normal.dotProduct(raydir))).add( color.scale(ks*std::max(float(0), normal.dotProduct(half))));
    return diffuse.add(specular);

    // return Vec3f(texture.red, texture.green, texture.blue);
}

void render(const std::vector<Triangle*> &triangle_list){

    int width = 200, height = 200;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);

    // read texture map
    png::image< png::rgb_pixel > texture_map("texture_map.png");

    // Trace rays
    int cnt = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -2);
            raydir.normalize();
            *pixel = trace(Vec3f(0,0.1,-7), raydir, triangle_list, texture_map);
        }
        std::cout << y << "\n"; cnt++;
    }
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./trial10.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
               (unsigned char)(std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
    delete [] image;
}


int main(){

    std::ifstream objinfile("input.obj");

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
