/*
 *  vlstructs.h
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/18/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *  This class contains various structures used in vrl,vpl,vsl,vbl rendering.
 */


#ifndef PBRT_VLSTRUCTS_H
#define PBRT_VLSTRUCTS_H
#include "pbrt.h"
#include "geometry.h"
#include "integrator.h"




// VirtualLight struct
struct VirtualLight {
    VirtualLight() { }
    VirtualLight(const Point &pp, const Normal &nn, const Spectrum &c,
                 float reps, BSDF *bs)
	: p(pp), n(nn), pathContrib(c), rayEpsilon(reps) { }
	
    Point p;
    Normal n;
    Spectrum pathContrib;
	
    float rayEpsilon;
	string toString(){
		std::ostringstream ss;
		ss<<"[";
		ss<<p.x;
		ss<<",";
		ss<<p.y;
		ss<<",";
		ss<<p.z;
		ss<<"]";
		return ss.str();	
	}
};

struct VirtualSphericalLight {
	VirtualSphericalLight(){}
	//missing incomming direction
	VirtualSphericalLight(const Point pp,const Vector incommingDirection, const Normal nn, const Spectrum c,
                 float reps,float rad, BSDF *bs)
	: p(pp),i(incommingDirection), n(nn), pathContrib(c), rayEpsilon(reps),radius(rad),bsdf(bs) { }
	BSDF *bsdf; // saving the lights bsdf  
	//added toString
	string toString(){
		std::ostringstream ss;
		ss<<"[";
		ss<<p.x;
		ss<<",";
		ss<<p.y;
		ss<<",";
		ss<<p.z;
		ss<<"]";
		return ss.str();	
	}
    string toPointString(){
		std::ostringstream ss;
		ss<<p.x;
		ss<<",";
		ss<<p.y;
		ss<<",";
		ss<<p.z;
		return ss.str();
	}
	//aka vl
	Point p;
    Normal n;
    //chybi mi tu incoming direction?
    Vector i;
    Spectrum pathContrib;
	
    float rayEpsilon;
	//MC added radius
	float radius;
};

inline string toPointString(VirtualSphericalLight &vsl){
    std::ostringstream ss;
    ss<<vsl.p.x;
    ss<<",";
    ss<<vsl.p.y;
    ss<<",";
    ss<<vsl.p.z;
    return ss.str();
}


struct VolumeVSL {
	VolumeVSL();
	
};

struct VolumePath{
    VolumePath(RayDifferential r){
        ray=r;
    }
    RayDifferential ray;
    //points with corresponding attenuated contribs
    std::vector<float> dists;
    //contribution of light
    Spectrum contrib;
    float getTransmittance(float distFromRayO){
        if (distFromRayO<0) {
            return 0;
        }
        //compute how many samples go beyond our sample point
        for (int i=0; i<dists.size(); ++i) {
            if (dists[i]>distFromRayO) {
                return (float)(dists.size()-(i+1))/(float)dists.size();
            }
        }
        return 0.;
    }
};

inline bool toLineString(VolumePath &vpth,string &data){
    if (vpth.dists.size()==0) {
        return false;
    }
    std::ostringstream ss;
    
    // Vector d=Normalize(vpth.ray.d);
    float lastDist=0.;
    Vector d=vpth.ray.d;
    Point p=vpth.ray.o;
    float distNorm=1./vpth.ray.maxt;
    //printf("\n=======\n");
    ss<<p.x;
    ss<<",";
    ss<<p.y;
    ss<<",";
    ss<<p.z;
    ss<<",";
    ss<<lastDist;
    ss<<",";
    for (int i=0; i<vpth.dists.size()-1; i++) {
        p=vpth.ray.o +d*vpth.dists[i];
        lastDist=vpth.dists[i]*distNorm;
        //printf("%f,",vpth.dists[i]);
        ss<<p.x;
        ss<<",";
        ss<<p.y;
        ss<<",";
        ss<<p.z;
        ss<<",";
        ss<<lastDist;
        ss<<",";
    }
    if (vpth.dists.size()>0) {
        int i=vpth.dists.size()-1;
        p=vpth.ray.o +d*vpth.dists[i];
        lastDist=vpth.dists[i]*distNorm;
        ss<<p.x;
        ss<<",";
        ss<<p.y;
        ss<<",";
        ss<<p.z;
        ss<<",";
        ss<<lastDist;
        ss<<",";
    }
    //and add the end point
    p=vpth.ray.o + d*vpth.ray.maxt;
    ss<<p.x;
    ss<<",";
    ss<<p.y;
    ss<<",";
    ss<<p.z;
    ss<<",";
    ss<<lastDist;
    data= ss.str();
    return true;
}

struct SurfaceLight{};


#endif