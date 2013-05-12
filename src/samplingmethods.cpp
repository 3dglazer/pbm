//
//  samplingmethods.cpp
//  pbrt
//
//  Created by Zdenek Glazer on 12.05.13.
//
//

#include "samplingmethods.h"


float A(float x, float phi, float h){
    return asinhf((x/h)*sinf(phi));
}

// analytical marginal inverse CDF for inverse squared sampling points vrl, phi is angle between camera vector and vrl vector
float invSqrCDF(float r,float phi,float h, float v0, float v1){
    return h*sinhf(Lerp(r,A(v0, phi, h), A(v1,phi,h)))/sinf(phi);
}

float B(float x,float h){
    return atanf(x/h);
}

//generates equiangular samples on the camera ray between u0 and u1 parametric distances, r is random number
float invEqAngCDF(float r,float h,float u0,float u1){
    return h*tanf(Lerp(r,B(u0,h), B(u1,h)));
}

//both rays have to be normalized + the vrl.mint is beginning of the media and vrl.maxt is the end of interraction. h will contain the smallest distance between the vectors
float sampleVRL(const Ray &r,const Ray &vrl, RNG& rng, float &h){
    float uh; //closest point parameter on camera ray
    float vh; //closest point parameter on vrl
    h=ln2lnPP(r, vrl, uh, vh);
    float phi=acosf(Dot(r.d, vrl.d));
    float v0d= vrl.mint-vh;
    float v1d= vrl.maxt-vh;
    float rnd=rng.RandomFloat();
    return invSqrCDF(rnd, phi, h, v0d, v1d);
}

//unit vectors expected otherwise wont work
float myline2line(const Point &p0,const Vector &v0,const Point &p1,const Vector &v1){
    bool paralel=false;
    Vector tmp=(v0-v1);
    if(tmp.Length()<0.00001) paralel=true;
    if (paralel) {
        tmp=Cross(Vector(p1-p0), v0);
        return tmp.Length();
    }else{
        tmp=Cross(v0,v1);
        return AbsDot(Vector(p1-p0), tmp)/tmp.Length();
    }
}

// Variation of the Code from http://geomalgorithms.com/a07-_distance.html
// References David Eberly, "Distance Methods" in 3D Game Engine Design (2006)
//Seth Teller, line_line_closest_points3d() (2000) cited in the Graphics  Algorithms FAQ (2001)
// Returns closest distance between lines and creates two closest points on these lines.
float ln2lnPP(const Ray &L1,const Ray &L2, float &p1, float &p2)
{
    Vector u = L1.d;
    Vector v = L2.d;
    Vector w = L1.o - L2.o;
    float a = Dot(u,u);
    float b = Dot(u,v);
    float c = Dot(v,v);
    float d = Dot(u,w);
    float e = Dot(v,w);
    float D = a*c - b*b;
    float sc, tc;
    if (D < 0.000001) {          // parallel
        sc = 0.0;
        tc = (b>c ? d/b : e/c);
    }
    else {
        sc = (b*e - c*d) / D;
        tc = (a*e - b*d) / D;
    }
    p1=sc;
    p2=tc;
    Vector ln=w+(sc*u) -(tc*v);
    return ln.Length();
}

//float
//dist3D_Line_to_Line( Line L1, Line L2)
//{
//    Vector   u = L1.P1 - L1.P0;
//    Vector   v = L2.P1 - L2.P0;
//    Vector   w = L1.P0 - L2.P0;
//    float    a = dot(u,u);         // always >= 0
//    float    b = dot(u,v);
//    float    c = dot(v,v);         // always >= 0
//    float    d = dot(u,w);
//    float    e = dot(v,w);
//    float    D = a*c - b*b;        // always >= 0
//    float    sc, tc;
//    
//    // compute the line parameters of the two closest points
//    if (D < SMALL_NUM) {          // the lines are almost parallel
//        sc = 0.0;
//        tc = (b>c ? d/b : e/c);    // use the largest denominator
//    }
//    else {
//        sc = (b*e - c*d) / D;
//        tc = (a*e - b*d) / D;
//    }
//    
//    // get the difference of the two closest points
//    Vector   dP = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)
//    
//    return norm(dP);   // return the closest distance
//}


float p2Ray(const Point &p,const Ray &v, float &vParam){
    Vector w = p - v.o;
    
    float c1 = Dot(w,v.d);
    float c2 = Dot(v.d,v.d);
    float b = c1 / c2;
    vParam=b;
    Point pnew = v.o + b * v.d;
    return Distance(p, pnew);
}

float p2Ray(const Point &p,const Point &vStart, const Point &vEnd,float &vParam){
    Vector v=Vector(vEnd-vStart);
    Vector w = p - vStart;
    
    float c1 = Dot(w,v);
    float c2 = Dot(v,v);
    float b = c1 / c2;
    vParam=b;
    Point pnew = vStart + b * v;
    return Distance(p, pnew);
}

float point2Segm( const Point &p,const Point &s0, const Point &s1, float &sParam){
    Vector v = s1 - s0;
    Vector w = p - s0;
    float c1 = Dot(w,v);
    if ( c1 <= 0 )
        sParam=0.;
    return Distance(p, s0);
    float c2 = Dot(v,v);
    if ( c2 <= c1 )
        sParam=1.;
    return Distance(p, s1);
    float b = c1 / c2;
    Point pn = s0 + b * v;
    sParam=b;
    return Distance(p, pn);
}

//===================================================================
