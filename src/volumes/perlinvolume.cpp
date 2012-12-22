/*
 *  perlinvolume.cpp
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/14/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */
#include "stdafx.h"
#include "volumes/perlinvolume.h"
#include "paramset.h"


// ExponentialDensity Method Definitions
PerlinDensity *CreatePerlinVolumeRegion(const Transform &volume2world,
												  const ParamSet &params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.);
    Spectrum Le = params.FindOneSpectrum("Le", 0.);
	
	//need number of octaves
	int nOctaves=params.FindOneInt("nOctaves", 3);
	float omega=params.FindOneFloat("omega", 0.5);
	float frequency=params.FindOneFloat("frequency", 2.0);
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
	bool inverted=params.FindOneBool("inverted", false);
    return new PerlinDensity(sigma_a, sigma_s, g, Le, BBox(p0, p1),
								  volume2world,nOctaves,omega,frequency,inverted);
}


