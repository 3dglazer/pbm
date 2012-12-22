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

class ParticleShooter {
public:	
	ParticleShooter(int rngSeed=1){seed=rngSeed;};
	~ParticleShooter(){
		volumePaths.empty();
		vsls.empty();
		psArena.FreeAll();
	};
	void shootParticles(const Scene * scene, Camera * camera, const Renderer *renderer, int nPaths,float radius=0.01);
	vector<VirtualSphericalLight> vsls;
private:
	int seed;
	vector<VolumePath> volumePaths;
	MemoryArena psArena;
};

#endif