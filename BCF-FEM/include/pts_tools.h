#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <random>
#include "Common.h"
//#include "ImageIO.h"
#include "Matrix3D.h"

#include <cmath>
#include <string.h>
#include <float.h>




//typedef struct Pt3D {
//    float x, y, z, w;
//} Pt3D;

//typedef struct Element {
//    int nodeIDs[3];
//} Element;

//typedef struct extrinsicParms {
//    float tx, ty, tz;
//    float thetaX, thetaY, thetaZ;
//} extrinsicParms;

void rotateLiver();
void generateIntraopSurface();
void addGaussianNoise(Pt3D* ptCloud, int nPts, float sigmaX, float sigmaY, float sigmaZ);
void rigidMotionTransformation(Pt3D* ptCloud, int nPts, extrinsicParms parms);
float ScTP(const Pt3D& a, const Pt3D& b, const Pt3D& c);
float VolTet(const Pt3D& a, const Pt3D& b, const Pt3D& c);
Pt3D crossProduct(const Pt3D& b, const Pt3D& c);
Pt3D Subtract(const Pt3D& a, const Pt3D& b);
void ConvertXYZtoOFF();