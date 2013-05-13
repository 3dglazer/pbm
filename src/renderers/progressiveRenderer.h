
/*
 pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.
 
 This file is part of pbrt.
 
 pbrt is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.  Note that the text contents of
 the book "Physically Based Rendering" are *not* licensed under the
 GNU GPL.
 
 pbrt is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 extended by Zdenek Glazer
 
 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_PROGRESSIVERENDERER_H
#define PBRT_RENDERERS_PROGRESSIVERENDERER_H

// renderers/samplerrenderer.h*
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
#include "particleshooter.h"
//MC BRUTEFORCE ON OFF
//#define BRUTE

// PGRenderer Declarations
class PGRenderer:public ProgressiveRenderer {
public:
    // SamplerRenderer Public Methods, should add something like nIterations,
	//MC added nIterations
    PGRenderer(Sampler *s, Camera *c,Camera *prc,Film* surface2surface,Film* surface2media,Film* media2surface, Film* media2media, ProgressiveSurfaceIntegrator *si,
                    ProgressiveVolumeIntegrator *vi, bool visIds,int nIterations,int nps,float rad,int rndSeed);
    ~PGRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
                const Sample *sample, RNG &rng, MemoryArena &arena,
                Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Lms(const Scene *scene, const RayDifferential &ray, const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect, Spectrum *T) const ;
    Spectrum Lss(const Scene *scene,
                 const RayDifferential &ray, const Sample *sample, RNG &rng,
                 MemoryArena &arena, Intersection *isect, Spectrum *T) const ;
    Spectrum Lsm(const Scene *scene,
                 const RayDifferential &ray, const Sample *sample, RNG &rng,
                 MemoryArena &arena, Intersection *isect, Spectrum *T) const ;
    Spectrum Lmm(const Scene *scene,
                 const RayDifferential &ray, const Sample *sample, RNG &rng,
                 MemoryArena &arena, Intersection *isect, Spectrum *T) const ;
    Spectrum Transmittance(const Scene *scene,
                           const RayDifferential &ray, const Sample *sample,
                           RNG &rng, MemoryArena &arena) const;
    float freeFlight(const Scene *scene, const Ray &r,Spectrum& tau,const RNG &rng) const;
private:
    // PGRenderer Private Data
    bool visualizeObjectIds;
    Sampler *sampler;
    Camera *camera;
    ProgressiveSurfaceIntegrator *surfaceIntegrator;
    ProgressiveVolumeIntegrator *volumeIntegrator;
	//MC added nIterations for rendering
	void renderIter(int currentIter,const Scene *scene,Sample* sample);
	int nIter;
	int seed;
	//MC added these params
	Camera* progCamera;
	//custom films for different light contribution computation types
	// surface to surface
	Film* ss;
	// media to surface
	Film* ms;
	// media to media
	Film* mm;
	// surface to media
	Film* sm;
	//MC particle shooter stuff 
	//ParticleShooter *particleShooter;
	int nParticles;
	float radius;
	//end of MC
};



// PGRendererTask Declarations
class PGRendererTask : public Task {
public:
    // PGRendererTask Public Methods added prc Camera which is used to compute averaged image
    PGRendererTask(const Scene *sc, ProgressiveRenderer *ren,ProgressiveVolumeIntegrator *vlmint, ProgressiveSurfaceIntegrator *surfint, Camera *c, Camera * prc,Film* surface2surface,Film* surface2media,Film* media2surface, Film* media2media,
                        ProgressReporter &pr, Sampler *ms, Sample *sam, 
                        bool visIds, int tn, int tc,int sd)
	: reporter(pr)
    {
        this->surfaceIntegrator=surfint;
        this->volumeIntegrator=vlmint;
        this->scene = sc; renderer = ren; camera = c; progCamera=prc; mainSampler = ms;
        this->origSample = sam; visualizeObjectIds = visIds; taskNum = tn; taskCount = tc;
		this->ss = surface2surface;
		this->ms = media2surface;
		this->mm = media2media;
		this->sm = surface2media;
		this->seed = sd;
		
    }
    void Run();
private:
    void repairSample(Spectrum &s);
    // PGRendererTask Private Data
    const Scene *scene;
    const ProgressiveRenderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    ProgressiveSurfaceIntegrator *surfaceIntegrator;
    ProgressiveVolumeIntegrator *volumeIntegrator;
    Sample *origSample;
    bool visualizeObjectIds;
    int taskNum, taskCount;
	//MC added these params
	Camera* progCamera;
	//custom films for different light contribution computation types
	// surface to surface
	Film* ss;
	// media to surface
	Film* ms;
	// media to media
	Film* mm;
	// surface to media
	Film* sm;
	
	int seed;
};



#endif // PBRT_RENDERERS_PGRenderer_H
