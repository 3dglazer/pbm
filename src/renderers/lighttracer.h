/*
 *  lighttracer.h
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/28/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_LIGHTTRACERRENDERER_H
#define PBRT_RENDERERS_LIGHTTRACERRENDERER_H
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
// will have to create lt task to go multithreaded
class LightTracerRenderer : public Renderer {
public:
	LightTracerRenderer(Camera *c,int nIterations,int pathsPerIter,int rndseed, SurfaceIntegrator *si, VolumeIntegrator *vi);
	~LightTracerRenderer();
	void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
				const Sample *sample, RNG &rng, MemoryArena &arena,
				Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
						   const Sample *sample, RNG &rng, MemoryArena &arena) const;
    #ifdef FREEFLIGHTEXTENSION
    float freeFlight(const Scene *scene, const Ray &r,Spectrum& tau,const RNG &rng) const;
#endif
	Camera* camera;
	int seed;
	int niter;
	float time;
	uint32_t maxPathCount;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
};

class LightShootingTask : public Task {
public:
	//musim se kouknout co vsechno potrebuju, zrejme scene a kamera
	LightShootingTask(const Scene* sscene, Camera* cm,Distribution1D *distrib,float tm,uint32_t maxPaths,int rndseed);
	void Run();
	// this method projects 3D point to camera and if it falls in film plate adds spectrum responce to the film.
	void filmAddSample(Point &p, Spectrum &sp);
	const Scene* scene;
	Camera* camera;
	int maxDepth;
	int seed;
	Distribution1D *lightDistribution;
	float time;
	uint32_t maxPathCount;
	Matrix4x4* world2Camera;
};
#endif