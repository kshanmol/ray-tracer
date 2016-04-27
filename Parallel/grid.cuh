#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#define HD __host__ __device__

struct Voxel{

    HD unsigned size() const { return voxelListSize; }
    HD Voxel() { ptr = 0; }
    HD Voxel(Triangle* op) {

		ptr = 0;
        triangleList = (Triangle*)malloc(sizeof(Triangle));
		triangleList = op;
		ptr++;

    }

    HD void AddPrimitive(Triangle* triangle) {
        triangleList[ptr] = *triangle;
		ptr++;
    }
    HD bool Intersect(const Ray &ray, Intersection *isect, Triangle& triangle_near, float& t0);


public:
	Triangle* triangleList;
	int voxelListSize;
	int ptr;

};

HD bool Voxel::Intersect(const Ray &ray, Intersection *isect, Triangle& triangle_near, float& t_min){

    bool hitSomething = false;

    for (int i = 0; i < voxelListSize; ++i) {
        Triangle* prim = &triangleList[i];

        if (prim->Intersect(ray, isect, triangle_near, t_min)){
			hitSomething = true;
        }
    }
    return hitSomething;

}

typedef struct Voxel Voxel;

class GridAccel{
public:
	//Grid accel public data
    Triangle* triangleList;
	int tListSize;
	int nVoxels[3];
	boundingBox bounds;
	float width[3], invWidth[3];
	Voxel **voxels;

	GridAccel(Triangle* triangle_list, int listSize);

    HD boundingBox WorldBound() const;
    HD bool CanIntersect() const { return true; }
    HD ~GridAccel();


	HD bool Intersect(const Ray &ray, Intersection *isect, Triangle& triangle_near, float& t0, int isDebugThread) const;

	HD bool IntersectP(const Ray &ray) const;	

private:

	HD int posToVoxel(const Vec3f& pt, int axis) const {

		float P[3] = {pt.x, pt.y, pt.z};
		float LB[3] = {bounds.lowerB.x, bounds.lowerB.y, bounds.lowerB.z};
		int v = (int)((P[axis] - LB[axis])*invWidth[axis]);
		return clamp(v, 0, nVoxels[axis] - 1);

	}

	HD float voxelToPos(int p, int axis) const {
		float lowerBound[3] = {bounds.lowerB.x, bounds.lowerB.y, bounds.lowerB.z};
		return lowerBound[axis] + p * width[axis];
	}

	HD inline int offset(int x, int y, int z) const {
		return z*nVoxels[0]*nVoxels[1] + y*nVoxels[0] + x;
	}

};

GridAccel::GridAccel(Triangle* triangle_list, int listSize){

	triangleList = (Triangle*)malloc(sizeof(Triangle)*tListSize);
	triangleList = triangle_list;
	tListSize = listSize;

	for(int i=0;i<tListSize;i++){

		bounds.union_(triangleList[i].v0);
		bounds.union_(triangleList[i].v1);
		bounds.union_(triangleList[i].v2);
	}

	Vec3f diff = bounds.upperB.subtract(bounds.lowerB);
	float delta[3] = {diff.x, diff.y, diff.z};

	int maxAxis = bounds.maxAxis();
	float maxInvWidth = 1.f/delta[maxAxis];

    float cubeRoot = 3.f * powf(float(tListSize), 1.f/3.f);

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

	int voxelTCount[totalVoxels];

	memset(voxelTCount, 0, totalVoxels* sizeof(int));

    //add primitives to grid voxels

    for(int i = 0;i<tListSize;i++){
    	//Find voxel extent of primitive

    	boundingBox* pb = new boundingBox();
    	pb->union_(triangleList[i].v0);
       	pb->union_(triangleList[i].v1);
    	pb->union_(triangleList[i].v2);

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
					voxelTCount[index]++;

    			}
    		}
    	}

    }

	for(int i = 0;i<totalVoxels;i++){
	
		if(voxelTCount[i] != 0){

			voxels[i] = new Voxel();
			voxels[i]->voxelListSize = voxelTCount[i];
			voxels[i]->triangleList = (Triangle*)malloc(voxelTCount[i] * sizeof(Triangle));
		}
	}

    for(int i = 0;i<tListSize;i++){
        //Find voxel extent of primitive

        boundingBox* pb = new boundingBox();
        pb->union_(triangleList[i].v0);
        pb->union_(triangleList[i].v1);
        pb->union_(triangleList[i].v2);

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
                 	
					if(voxelTCount[index] != 0)
                    	voxels[index]->AddPrimitive(&triangleList[i]);

                }
            }
        }

    }

    //Construction over
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

//From pbr
HD bool GridAccel::Intersect(const Ray& ray, Intersection *isect, Triangle& triangle_near, float& t_min, int isDebugThread) const{

	//Check ray against overall grid bounds
	float rayT;

	if(bounds.Inside(ray(ray.mint)))
	 	rayT = ray.mint;
	else if (!bounds.Intersect(ray, isDebugThread, &rayT)){
		return false;
	}

	Vec3f gridIntersect = ray(rayT);
	float grid_intersect[3] = {gridIntersect.x, gridIntersect.y, gridIntersect.z};

	//Setup walk through algorithm

    float NextCrossingT[3], DeltaT[3]; //Parametric position of next voxel, parametric width of ray along each axis
    int Step[3], Out[3], Pos[3]; //Stepping direction ( +1 or -1 ), End of grid marker, current voxel coordinates  
    float ray_dir[3] = {ray.raydir.x, ray.raydir.y, ray.raydir.z};

    for (int axis = 0; axis < 3; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray_dir[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT + (voxelToPos(Pos[axis]+1, axis) - grid_intersect[axis]) / ray_dir[axis];
            DeltaT[axis] = width[axis] / ray_dir[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT + (voxelToPos(Pos[axis], axis) - grid_intersect[axis]) / ray_dir[axis];
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
       
		if (voxel != NULL){
            hitSomething |= voxel->Intersect(ray, isect, triangle_near, t_min);
		}
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

