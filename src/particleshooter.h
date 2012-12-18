/*
 *  particleshooter.h
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/18/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */
#ifndef PBRT_PARTICLESHOOTER_H
#define PBRT_PARTICLESHOOTER_H
#include "pbrt.h"
#include "vlstructs.h"

class ParticleShooter {
public:	
	ParticleShooter(int rngSeed=1){seed=rngSeed;};
	~ParticleShooter();
	shootParticles(const Scene *, const Camera *, const Renderer *,const int nParticles,const float radius=0.01);
private:
	int seed;
	vector<VolumePaths> volumePaths;
	vector<SurfaceLights> surfaceLights;
	MemoryArena psArena;
};

#endif