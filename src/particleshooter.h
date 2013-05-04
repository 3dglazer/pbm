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
//#include "vlstructs.h"
#include "datadumper.h"
#include "scene.h"
#include "camera.h"
#include "volume.h"
#include "rng.h"
#include <vector.h>

class ParticleShooter {
public:	
	ParticleShooter(int rngSeed=1,int maxScatteringEvents = 8){seed=rngSeed;maxScattering=maxScatteringEvents;};
	~ParticleShooter(){
		volumePaths.erase(volumePaths.begin(),volumePaths.end()); //ensure the objects will be deleted
		vsls.erase(vsls.begin(),vsls.end());
		psArena.FreeAll();
	};
	void shootParticles(const Scene * scene, Camera * camera, const Renderer *renderer, int nPaths,float radius=0.01);
	vector<VirtualSphericalLight *> vsls;
    std::vector<VolumePath *> volumePaths;
private:
	int seed;
    int maxScattering;
	MemoryArena psArena;
};

#endif