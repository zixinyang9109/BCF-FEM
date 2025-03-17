#include <fstream>
#include <cfloat>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Common.h"
using namespace std;

//https://github.com/ucanbizon/downsampling-point-cloud/blob/master/downsample.cpp//

class point_and_box
{
public:
    int idx;
    int box;
    point_and_box(int arg_idx) { idx = arg_idx; box = -1; }
    bool operator < (const point_and_box& rhs) const { return(box < rhs.box); }
};

template <class T> class vec3
{
public:
    T x, y, z;

    vec3() {}

    vec3(T arg_x, T arg_y, T arg_z)
    {
        x = arg_x;
        y = arg_y;
        z = arg_z;
    }

    vec3<T>& operator = (const vec3<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    vec3(const vec3<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
    }

    vec3<T>& operator+=(const vec3<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
};

void getMinMax(vector< vec3<double> >& inCloud, vec3<double>& minp, vec3<double>& maxp);


void filterPoints(vec3<float>& leafCount, vector< vec3<double> >& inCloud, vector< vec3<double> >& outCloud, vector< point_and_box > indices);


void downsamplePointCloud(Pt3D* ptCloud, int& nPts, float nleafCount, bool bFlag);
