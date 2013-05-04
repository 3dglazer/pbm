/*
 *  photonDisc.cpp
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/22/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */

#include "photonDisc.h"


PhotonDisc::PhotonDisc(Point* position,float rad){
	radius=rad;
	radiusSQ=rad*rad;
	pos=position;
}

float PhotonDisc::Area() const {
    return piFT*radius*radius*radius;
}

BBox PhotonDisc::ObjectBound() const {
    return BBox(Point(-radius, -radius, -radius),
                Point( radius,  radius, radius));
}

bool PhotonDisc::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                       DifferentialGeometry *dg) const {
	//r.d is ray direction thus a disc normal
	Vector rop=Vector(pos->x-r.o.x,pos->y-r.o.y,pos->z-r.o.z);
	float a=Dot(rop,r.d);
	float bSQ= a*a-rop.LengthSquared();
	if (bSQ>radiusSQ) {
		return false;
	}
	*tHit=a;
	// Compute _rayEpsilon_ for quadric intersection
    *rayEpsilon = 5e-4f * *tHit;
    return true;
}

bool PhotonDisc::IntersectP(const Ray &r) const {
	//r.d is ray direction thus a disc normal
	Vector rop=Vector(pos->x-r.o.x,pos->y-r.o.y,pos->z-r.o.z);
	float a=Dot(rop,r.d);
	float bSQ= a*a-rop.LengthSquared();
	if (bSQ<radiusSQ) {
		return true;
	}else {
		return false;
	}
}
