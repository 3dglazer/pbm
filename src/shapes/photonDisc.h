/*
 *  photonDisc.h
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/22/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_PHOTONDISC_H
#define PBRT_SHAPES_PHOTONDISC_H

#include "shape.h"

//MC chybi dedicnost od : public Shape
class PhotonDisc  {
public:
    // PhotonDisc Public Methods
	PhotonDisc();
    PhotonDisc(Point* position,float rad);
    BBox ObjectBound() const;
	float Area() const;
    bool IntersectP(const Ray &ray) const;
	bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
	
private:
    // Photon Disc
	Point* pos;
    float radius;
	float radiusSQ;
	//piFT=4/3*pi
	static const float piFT= 4.1887902047863905;
};


#endif 