
double det(double a1, double a2, double a3,
            double b1, double b2, double b3,
            double c1, double c2, double c3);

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

struct Intersection;

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


class boundingBox{

public:
    //public data
    Vec3f lowerB, upperB; 

    //public methods
    void union_(const Vec3f& vert);
    int maxAxis();
    boundingBox();
    bool Inside(const Vec3f& pt) const;
    bool Intersect(const Ray &ray, int isDebugThread, float *hitt0 = NULL, float *hitt1 = NULL) const;

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

inline void debug(Vec3f v){
    printf("%f %f %f\n", v.x, v.y, v.z);
}

bool boundingBox::Intersect(const Ray &ray, int isDebugThread, float *hitt0 , float *hitt1 ) const {

    float t0 = ray.mint, t1 = ray.maxt;

    float minB[3] = {lowerB.x, lowerB.y, lowerB.z};
    float maxB[3] = {upperB.x, upperB.y, upperB.z};
    float ray_orig[3] = {ray.orig.x, ray.orig.y, ray.orig.z};
    float ray_dir[3] = {ray.raydir.x, ray.raydir.y, ray.raydir.z};

    if(isDebugThread){
        debug(lowerB);
        debug(upperB);
        debug(ray.orig);
        debug(ray.raydir);
    }

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

