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
double det(double a1, double a2, double a3,
            double b1, double b2, double b3,
            double c1, double c2, double c3);
int clamp(int what, int low, int high);
class Triangle;

float global_t = INFINITY;
Triangle* global_triangle_near = NULL;

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

class Ray{

public:

	//Methods
    Ray() : mint(0.f), maxt(INFINITY), depth(0) { }
    Ray(const Vec3f &origin, const Vec3f &direction, float start, float end = INFINITY, int d = 0)
        : orig(origin), raydir(direction), mint(start), maxt(end), depth(d) { }
    Ray(const Vec3f &origin, const Vec3f &direction, const Ray &parent, float start, float end = INFINITY)
        : orig(origin), raydir(direction), mint(start), maxt(end), depth(parent.depth+1) { }

	Vec3f operator()(float t) const { return orig.add(raydir.scale(t)); }

	//Data
	Vec3f orig;
	Vec3f raydir;
	mutable float maxt, mint;
	int depth;

};

struct Intersection{



};

typedef struct Intersection Intersection;

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

    bool Intersect(const Ray& ray, Intersection *isect) const{

    	Vec3f orig = ray.orig, dir = ray.raydir;

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
            if(t_ < global_t){
            	global_t = t_;
            	global_triangle_near = (Triangle*)this;
            }
            return true;
        }
        return false;

    }

    bool rayTriangleIntersect(Vec3f &orig, const Vec3f &dir, float &t, float &beta, float &gamma){

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

struct Voxel{

    unsigned size() const { return triangleList.size(); }
    Voxel() { }
    Voxel(Triangle* op) {
        triangleList.push_back(op);
    }

    void AddPrimitive(Triangle* triangle) {
        triangleList.push_back(triangle);
    }
    bool Intersect(const Ray &ray, Intersection *isect);


private:
	std::vector<Triangle*> triangleList;

};

bool Voxel::Intersect(const Ray &ray, Intersection *isect){

    bool hitSomething = false;
    for (int i = 0; i < triangleList.size(); ++i) {
        Triangle* prim = triangleList[i];
        if (prim->Intersect(ray, isect)){
            hitSomething = true;
        }
    }
    return hitSomething;

}

typedef struct Voxel Voxel;


class boundingBox{

public:
	//public data
	Vec3f lowerB, upperB; 

	//public methods
	void union_(const Vec3f& vert);
	int maxAxis();
	boundingBox();
	bool Inside(const Vec3f& pt) const;
	bool Intersect(const Ray &ray, float *hitt0 = NULL, float *hitt1 = NULL) const;

};

boundingBox::boundingBox(void){
	lowerB = Vec3f(INFINITY, INFINITY, INFINITY);
	upperB = Vec3f(-INFINITY, -INFINITY, -INFINITY);
}

void boundingBox::union_(const Vec3f& vert){

	upperB.x = std::max(vert.x, upperB.x);
	upperB.y = std::max(vert.y, upperB.y);
	upperB.z = std::max(vert.z, upperB.z);

	lowerB.x = std::min(vert.x, lowerB.x);
	lowerB.y = std::min(vert.y, lowerB.y);
	lowerB.z = std::min(vert.z, lowerB.z);

}

int boundingBox::maxAxis(){

	int axis = ( upperB.x - lowerB.x > upperB.y - lowerB.y ) ? 0 : 1;
	if(axis)
		axis = ( upperB.y - lowerB.y > upperB.z - lowerB.z ) ? 1 : 2;
	else
		axis = ( upperB.x - lowerB.x > upperB.z - lowerB.z ) ? 0 : 2;

	return axis;
}

bool boundingBox::Inside(const Vec3f &pt) const {
    return (pt.x >= lowerB.x && pt.x <= upperB.x && pt.y >= lowerB.y && pt.y <= upperB.y && pt.z >= lowerB.z && pt.z <= upperB.z);
}

bool boundingBox::Intersect(const Ray &ray, float *hitt0 , float *hitt1 ) const {

    float t0 = ray.mint, t1 = ray.maxt;

    float minB[3] = {lowerB.x, lowerB.y, lowerB.z};
    float maxB[3] = {upperB.x, upperB.y, upperB.z};
    float ray_orig[3] = {ray.orig.x, ray.orig.y, ray.orig.z};
    float ray_dir[3] = {ray.raydir.x, ray.raydir.y, ray.raydir.z};

    for (int i = 0; i < 3; ++i) {
        // Update interval for _i_th bounding box slab
        float invRayDir = 1.f / ray_dir[i];
        float tNear = (minB[i] - ray_orig[i]) * invRayDir;
        float tFar  = (maxB[i] - ray_orig[i]) * invRayDir;

        // Update parametric interval from slab intersections
        if (tNear > tFar) std::swap(tNear, tFar);
        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar  < t1 ? tFar  : t1;
        if (t0 > t1) return false;
    }
    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;
    return true;

}

class GridAccel{
public:
	//Grid accel public data
    std::vector<Triangle*> triangleList;
	int nVoxels[3];
	boundingBox bounds;
	float width[3], invWidth[3];
	Voxel **voxels;

	GridAccel(std::vector<Triangle*>& triangle_list);

    boundingBox WorldBound() const;
    bool CanIntersect() const { return true; }
    ~GridAccel();


	bool Intersect(const Ray &ray, Intersection *isect) const;

	bool IntersectP(const Ray &ray) const;	

private:

	int posToVoxel(const Vec3f& pt, int axis) const {

		float P[3] = {pt.x, pt.y, pt.z};
		float LB[3] = {bounds.lowerB.x, bounds.lowerB.y, bounds.lowerB.z};
		int v = (int)((P[axis] - LB[axis])*invWidth[axis]);
		return clamp(v, 0, nVoxels[axis] - 1);

	}

	float voxelToPos(int p, int axis) const {
		float lowerBound[3] = {bounds.lowerB.x, bounds.lowerB.y, bounds.lowerB.z};
		return lowerBound[axis] + p * width[axis];
	}

	inline int offset(int x, int y, int z) const {
		return z*nVoxels[0]*nVoxels[1] + y*nVoxels[0] + x;
	}

};

GridAccel::GridAccel(std::vector<Triangle*>& triangle_list){

	triangleList = triangle_list;
	for(int i=0;i<triangleList.size();i++){
		bounds.union_(triangleList[i]->v0);
		bounds.union_(triangleList[i]->v1);
		bounds.union_(triangleList[i]->v2);
	}

	Vec3f diff = bounds.upperB.subtract(bounds.lowerB);
	float delta[3] = {diff.x, diff.y, diff.z};

	int maxAxis = bounds.maxAxis();
	float maxInvWidth = 1.f/delta[maxAxis];

    float cubeRoot = 3.f * powf(float(triangleList.size()), 1.f/3.f);

    float voxelsPerUnitDist = cubeRoot * maxInvWidth;

    for (int axis = 0; axis < 3; ++axis) {
        nVoxels[axis] = (int) (delta[axis] * voxelsPerUnitDist + 1);
        nVoxels[axis] = clamp(nVoxels[axis], 1, 64);
    }

	//Compute voxel widths
    for (int axis = 0; axis < 3; ++axis) {
        width[axis] = delta[axis] / nVoxels[axis];
        invWidth[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
    }

    //Allocate memory to voxels
    int totalVoxels = nVoxels[0] * nVoxels[1] * nVoxels[2];

    voxels = (struct Voxel* *)malloc(totalVoxels*sizeof(struct Voxel*));

    memset(voxels, 0, totalVoxels * sizeof(struct Voxel*));

    //add primitives to grid voxels

    for(int i = 0;i<triangleList.size();i++){
    	//Find voxel extent of primitive

    	boundingBox* pb = new boundingBox();
    	pb->union_(triangleList[i]->v0);
       	pb->union_(triangleList[i]->v1);
    	pb->union_(triangleList[i]->v2);

    	int vmin[3], vmax[3];

    	for(int axis = 0;axis < 3;axis++){
			vmin[axis] = posToVoxel(pb->lowerB, axis);
			vmax[axis] = posToVoxel(pb->upperB, axis);
    	}

    	//Add primitive to overlapping voxels

    	for(int x = vmin[0];x<=vmax[0];x++){
    		for(int y = vmin[1];y<=vmax[1];y++){
    			for(int z = vmin[2];z<=vmax[2];z++){

    				int index = offset(x,y,z);

    				if(!voxels[index])
    					voxels[index] = new Voxel(triangleList[i]);
    				else
    					voxels[index]->AddPrimitive(triangleList[i]);

    			}
    		}
    	}

    }

    //Construction over

    for(int i = 0;i<3;i++)
    	std::cout << nVoxels[i] << " ";
	std::cout << "\n";

	int count = 0;
	for(int i = 0;i<totalVoxels;i++){
		if(voxels[i])
			if(voxels[i]->size() > 5)
				count++;
	}
	std::cout << count << "\n";

}

boundingBox GridAccel::WorldBound() const {
    return bounds;
}

//destructor
GridAccel::~GridAccel() {
    for (int i = 0; i < nVoxels[0]*nVoxels[1]*nVoxels[2]; ++i)
        if (voxels[i]) voxels[i]->~Voxel();
    free(voxels);
}

//copied 
bool GridAccel::Intersect(const Ray& ray, Intersection *isect ) const{

	//Check ray against overall grid bounds
	float rayT;

	if(bounds.Inside(ray(ray.mint)))
		rayT = ray.mint;
	else if (!bounds.Intersect(ray, &rayT))
		return false;

	Vec3f gridIntersect = ray(rayT);
	float grid_intersect[3] = {gridIntersect.x, gridIntersect.y, gridIntersect.z};

	//Setup 3D DDA for ray

    float NextCrossingT[3], DeltaT[3];
    int Step[3], Out[3], Pos[3];
    float ray_dir[3] = {ray.raydir.x, ray.raydir.y, ray.raydir.z};

    for (int axis = 0; axis < 3; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray_dir[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis]+1, axis) - grid_intersect[axis]) / ray_dir[axis];
            DeltaT[axis] = width[axis] / ray_dir[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis], axis) - grid_intersect[axis]) / ray_dir[axis];
            DeltaT[axis] = -width[axis] / ray_dir[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

	//Walk ray through voxel grid

    bool hitSomething = false;
    for (;;) {
        // Check for intersection in current voxel and advance to next
        Voxel *voxel = voxels[offset(Pos[0], Pos[1], Pos[2])];
        if (voxel != NULL)
            hitSomething |= voxel->Intersect(ray, isect);

        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int bits = ((NextCrossingT[0] < NextCrossingT[1]) << 2) +
                   ((NextCrossingT[0] < NextCrossingT[2]) << 1) +
                   ((NextCrossingT[1] < NextCrossingT[2]));
        const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
        int stepAxis = cmpToAxis[bits];
        if (ray.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
    return hitSomething;

}


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

Vec3f fast_trace(Ray& ray, const std::vector<Triangle*> &triangle_list, GridAccel* newGridAccel)
{

	Intersection* isect;
	Vec3f rayorig = ray.orig, raydir = ray.raydir;
	global_t = INFINITY;

	bool hitSomething = newGridAccel->Intersect(ray, isect);    
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

void render(const std::vector<Triangle*> &triangle_list, GridAccel* newGridAccel){

    int width = 256, height = 256;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.0);

    // Trace rays
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, 2);
            raydir.normalize();
            Ray ray(Vec3f(0,0.1,-7), raydir, 0);
            *pixel = trace(ray, Vec3f(0,0.1,-7), raydir, triangle_list);
            //*pixel = fast_trace(ray, triangle_list, newGridAccel);
        }
        //std::cout << y << "\n";
    }
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./start_comp.ppm", std::ios::out | std::ios::binary);
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

    double x_max = -INFINITY, x_min = INFINITY, y_max = -INFINITY, y_min = INFINITY, z_max = -INFINITY, z_min = INFINITY;

    for(int i = 0;i<vertices.size();i++){
    	if(vertices[i]->x > x_max)
    		x_max = vertices[i]->x;
       	if(vertices[i]->y > y_max)
    		y_max = vertices[i]->y;
       	if(vertices[i]->z > z_max)
    		z_max = vertices[i]->z;
       	if(vertices[i]->x < x_min)
    		x_min = vertices[i]->x;
       	if(vertices[i]->y < y_min)
    		y_min = vertices[i]->y;
       	if(vertices[i]->z < z_min)
    		z_min = vertices[i]->z;
    }

    double num_objects = triangle_list.size();

    double x_extent = x_max - x_min, y_extent = y_max - y_min, z_extent = z_max - z_min;

    double max_extent = std::max(x_extent, std::max(y_extent, z_extent));

    double cube_root = 3.0f*powf(num_objects, 1.f/3.f);

    std::cout << x_max << " " << x_min << " " << x_extent << "\n";
    std::cout << y_max << " " << y_min << " " << y_extent << "\n";
    std::cout << z_max << " " << z_min << " " << z_extent << "\n";

    std::cout << max_extent << " " << cube_root << "\n";

    GridAccel* newGridAccel = new GridAccel(triangle_list);

    std::cout << newGridAccel->bounds.lowerB << "\n";
    std::cout << newGridAccel->bounds.upperB << "\n";


    render(triangle_list, newGridAccel);

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
