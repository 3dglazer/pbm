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

// VirtualLight struct
struct VirtualLight {
    VirtualLight() { }
    VirtualLight(const Point &pp, const Normal &nn, const Spectrum &c,
                 float reps)
	: p(pp), n(nn), pathContrib(c), rayEpsilon(reps) { }
	
    Point p;
    Normal n;
    Spectrum pathContrib;
	
    float rayEpsilon;

};

struct VirtualSphericalLight {
	VirtualSphericalLight(){}
	
	VirtualSphericalLight(const Point &pp, const Normal &nn, const Spectrum &c,
                 float reps,float rad, BSDF *bs)
	: p(pp), n(nn), pathContrib(c), rayEpsilon(reps),radius(rad),bsdf(bs) { }
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
	//aka vl
	Point p;
    Normal n;
    Spectrum pathContrib;
	
    float rayEpsilon;
	//MC added radius
	float radius;
}


#endif