/*
 *  perlinvolume.h
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/14/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_VOLUMES_PERLINVOLUME_H
#define PBRT_VOLUMES_PERLINVOLUME_H

#include "volume.h"
#include "texture.h"

// ExponentialDensity Declarations
class PerlinDensity : public DensityRegion {
public:
    // ExponentialDensity Public Methods
    PerlinDensity(const Spectrum &sa, const Spectrum &ss,
                       float gg, const Spectrum &emit, const BBox &e,
                       const Transform &v2w,int nOctaves,float omega)
	: DensityRegion(sa, ss, gg, emit, v2w), extent(e) ,omg(omega),octaves(nOctaves){
    }
    BBox WorldBound() const { return Inverse(WorldToVolume)(extent); }
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
    float Density(const Point &Pobj) const {
        if (!extent.Inside(Pobj)) return 0;
        return SimpleTurbulence(Pobj, omg, octaves);
    }
private:
    // PerlinDensity Private Data
    BBox extent;
	float omg;
	int octaves;
};


PerlinDensity *CreatePerlinVolumeRegion(const Transform &volume2world,
												  const ParamSet &params);

#endif // PBRT_VOLUMES_EXPONENTIAL_H
