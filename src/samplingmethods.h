//
//  samplingmethods.h
//  pbrt
//
//  Created by Zdenek Glazer on 12.05.13.
//
//

#ifndef __pbrt__samplingmethods__
#define __pbrt__samplingmethods__

#include <iostream>
//#include <cmath>
#include <math.h>
#include "pbrt.h"
#include <stdlib.h>
#include "primitive.h"
#include "rng.h"
float A(float x, float phi, float h);
float invSqrCDF(float r,float phi,float h, float v0, float v1);
float B(float x,float h);
float invEqAngCDF(float r,float h,float u0,float u1);
float sampleVRL(const Ray &r,const Ray &vrl, RNG& rng, float &h);
float myline2line(const Point &p0,const Vector &v0,const Point &p1,const Vector &v1);
float ln2lnPP(const Ray &L1,const Ray &L2, float &p1, float &p2);

//p the point, s0 segment start, s1 segment endPoint, the parametric distance to the point on segment will be returned in sParam;
float point2Segm( const Point &p,const Point &s0, const Point &s1, float &sParam);

float p2Ray(const Point &p,const Ray &v, float &vParam);
float p2Ray(const Point &p,const Point &vStart, const Point &vEnd,float &vParam);
#endif /* defined(__pbrt__samplingmethods__) */
