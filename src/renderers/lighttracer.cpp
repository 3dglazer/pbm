/*
 *  lighttracer.cpp
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/28/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */

#include "renderers/lighttracer.h"
#include "stdafx.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "montecarlo.h"
LightShootingTask::LightShootingTask(const Scene* sscene,Camera* cm,Distribution1D *distrib,float tm,uint32_t maxPaths,int rndseed){
	camera=cm;
	seed=rndseed;
	scene=sscene;
	lightDistribution=distrib;
	time=tm;
	maxPathCount=maxPaths;
	Transform* tr= new Transform();
	camera->CameraToWorld.Interpolate(time, tr);
	//for 3d point projection , its without motion blur
	world2Camera=new Matrix4x4(tr->GetInverseMatrix());
	printf("%s",world2Camera->toString().c_str());
	delete tr;
}

void LightShootingTask::filmAddSample(Point &p, Spectrum &sp){
	//projektnu bod 
	// asi by slo udelat i motion blur:
	//camera->CameraToWorld.Interpolate(time,transf)
}

void LightShootingTask::Run(){
	// tady by mel byt kod z photon mappingu
	
	MemoryArena arena;
	uint32_t totalPaths = 0;
    RNG rng(seed);
	 PermutedHalton halton(6, rng);
	while (true) {
        // Follow photon paths for a block of samples
        const uint32_t blockSize = 4096;
        for (uint32_t i = 0; i < blockSize; ++i) {
            float u[6];
            halton.Sample(++totalPaths, u);
			// Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];
			
            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],time, &photonRay, &Nl, &pdf);
            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
			
            if (!alpha.IsBlack()) {
                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);
                bool specularPath = true;
                Intersection photonIsect;
                int nIntersections = 0;
                while (scene->Intersect(photonRay, &photonIsect)) {
                    ++nIntersections;
					//MC tady by mel byt i kod pro volumetriku
					
                    // Handle photon/surface intersection
                   // alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
                    BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
					
                    Vector wo = -photonRay.d;
					//MC tady se ukladaly photony takze tady bych mel ukladat samples do filmu kamery
					//  // Deposit photon at surface
					//Photon photon(photonIsect.dg.p, alpha, wo);
					//tuhle metodu chci pouzit
					//filmAddSample()
					
					if (nIntersections >= maxDepth) break;
					
                    // Sample new photon ray direction
                    Vector wi;
                    float pdf;
                    BxDFType flags;
                    Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
                                                       &pdf, BSDF_ALL, &flags);
					
                    if (fr.IsBlack() || pdf == 0.f) break;
                    Spectrum anew = alpha * fr *
					AbsDot(wi, photonBSDF->dgShading.nn) / pdf;
					
                    // Possibly terminate photon path with Russian roulette
                    float continueProb = min(1.f, anew.y() / alpha.y());
                    if (rng.RandomFloat() > continueProb)
                        break;
                    alpha = anew / continueProb;
                    specularPath &= ((flags & BSDF_SPECULAR) != 0);
                    
                    photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
                                                photonIsect.rayEpsilon);
					
				}
                PBRT_PHOTON_MAP_FINISHED_RAY_PATH(&photonRay, &alpha);
			}
			
		arena.FreeAll();
		}
		
		//termination criteria ???
		if (totalPaths==maxPathCount) {
			break;
		}
	}
}

LightTracerRenderer::LightTracerRenderer(Camera *c,int nIterations,int pathsPerIter,int rndseed, SurfaceIntegrator *si, VolumeIntegrator *vi){
	camera=c;
	seed=rndseed;
	niter=nIterations;
	maxPathCount=pathsPerIter;
    
}
LightTracerRenderer::~LightTracerRenderer() {
	delete camera;
}

void LightTracerRenderer::Render(const Scene *scene) {
	// Compute light power CDF for photon shooting
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
	
	vector<Task *> lightShootingTasks;
    int nTasks = NumSystemCores();
	//iteratively increase fidelity of the image
	for (int iter=0; iter<niter; iter++) {
		//multi tasking one for each processor
		for (int i = 0; i < nTasks; ++i)
			lightShootingTasks.push_back(new LightShootingTask(scene,camera,lightDistribution, camera ? camera->shutterOpen : 0.f,maxPathCount,seed*(i+1)));
		EnqueueTasks(lightShootingTasks);
		WaitForAllTasks();
		for (uint32_t i = 0; i < lightShootingTasks.size(); ++i)
			delete lightShootingTasks[i];
		//write image
		camera->film->WriteImage();
	}

	
}

Spectrum LightTracerRenderer::Li(const Scene *scene,
							 const RayDifferential &ray, const Sample *sample, RNG &rng,
							 MemoryArena &arena, Intersection *isect, Spectrum *T) const {
	return NULL;
}

//MC
float LightTracerRenderer::freeFlight(const Scene *scene, const Ray &r,Spectrum& tau,const RNG &rng) const{
    return volumeIntegrator->freeFlight(scene, r, tau, rng);
}

Spectrum LightTracerRenderer::Transmittance(const Scene *scene,
										const RayDifferential &ray, const Sample *sample, RNG &rng,
										MemoryArena &arena) const {
	return NULL;

}