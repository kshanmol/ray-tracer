#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>

static const float eps = 1e-8; 

template<typename T>
class Vec3{
public:

    T x, y, z;
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

    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }

    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v){
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

typedef Vec3<float> Vec3f;

class Triangle{
public:
    
    Vec3f v0, v1, v2;

    Triangle(
        const Vec3f &v_0,
        const Vec3f &v_1, 
        const Vec3f &v_2 ) :
        v0(v_0), v1(v_1), v2(v_2)
    { /* empty */ }
    bool rayTriangleIntersect(const Vec3f &orig, const Vec3f &dir, float &t){

        Vec3f v0v1 = v1 - v0; 
        Vec3f v0v2 = v2 - v0; 
        Vec3f pvec = dir.crossProduct(v0v2); 
        float det = v0v1.dotProduct(pvec); 
     
        // ray and triangle are parallel if det is close to 0
        if (fabs(det) < eps) return false; 
     
        float invDet = 1 / det; 
     
        Vec3f tvec = orig - v0; 
        float u = tvec.dotProduct(pvec) * invDet; 
        if (u < 0 || u > 1) return false; 
     
        Vec3f qvec = tvec.crossProduct(v0v1); 
        float v = dir.dotProduct(qvec) * invDet; 
        if (v < 0 || u + v > 1) return false; 
     
        t = v0v2.dotProduct(qvec) * invDet; 
     
        return true; 
    }

};

Vec3f trace(const Vec3f &rayorig, const Vec3f &raydir, const std::vector<Triangle*> &triangle_list){

    float tnear = INFINITY;

    const Triangle* triangle_near = NULL;
    for (int i = 0; i < triangle_list.size(); ++i) {
        float t0 = INFINITY;
        if (triangle_list[i]->rayTriangleIntersect(rayorig, raydir, t0)) {
            if (t0 < tnear) {
                tnear = t0;
                triangle_near = triangle_list[i];
            }
        }
    }
    if (!triangle_near)
        return Vec3f(0);

    std::cout << tnear << " " << "Found\n";
    return Vec3f(fabs(0.0625*tnear));
    
}

void render(const std::vector<Triangle*> &triangle_list){

    int width = 200, height = 200;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace rays
    int cnt = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -2);
            raydir.normalize();
            *pixel = trace(Vec3f(0,0.1,-7), raydir, triangle_list);
            std::cout << cnt << "\n"; cnt++;
        }
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
        else if (type_.compare("f") == 0){

            iss >> fv1 >> fv2 >> fv3;
            std::stringstream ssfv1(fv1);
            std::stringstream ssfv2(fv2);
            std::stringstream ssfv3(fv3);

            int v1, v2, v3;
            ssfv1 >> v1;
            ssfv2 >> v2;
            ssfv3 >> v3;

            triangle = new Triangle(*vertices[v1-1], *vertices[v2-1], *vertices[v3-1]);
            triangle_list.push_back(triangle);
        }
    }

    render(triangle_list);

    return 0;
}
