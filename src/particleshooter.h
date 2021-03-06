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
	ParticleShooter(int rngSeed=1,int maxScatteringEvents = 16){seed=rngSeed;maxScattering=maxScatteringEvents;};
	~ParticleShooter(){
		volumePaths.erase(volumePaths.begin(),volumePaths.end()); //ensure the objects will be deleted
		vsls.erase(vsls.begin(),vsls.end());
		psArena.FreeAll();
	};
	void shootParticles(const Scene * scene, Camera * camera, const ProgressiveRenderer *renderer, int nPaths,float radius=0.01);
    void dumpVSLS(std::string fileName){
        DataDumper dd;
        dd.dumpVirtualSphericalLight(vsls,fileName.c_str());
    }
    void dumpVpths(std::string fileName){
        DataDumper dd;
        dd.dumpVolumePaths(volumePaths,fileName.c_str());
    }
	vector<VirtualSphericalLight *> vsls;
    std::vector<VolumePath *> volumePaths;
private:
	int seed;
    int maxScattering;
	MemoryArena psArena;
};

#endif