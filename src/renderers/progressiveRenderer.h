
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

// ProgressiveRenderer Declarations
class ProgressiveRenderer : public Renderer {
public:
    // SamplerRenderer Public Methods, should add something like nIterations,
	//MC added nIterations
    ProgressiveRenderer(Sampler *s, Camera *c,Camera *prc, SurfaceIntegrator *si,
                    VolumeIntegrator *vi, bool visIds,int nIterations,int rndSeed);
    ~ProgressiveRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
				const Sample *sample, RNG &rng, MemoryArena &arena,
				Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
						   const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // ProgressiveRenderer Private Data
    bool visualizeObjectIds;
    Sampler *sampler;
    Camera *camera;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
	//MC added nIterations for rendering
	void renderIter(int currentIter,const Scene *scene,Sample* sample);
	int nIter;
	int seed;
	Camera* progCamera;
};



// ProgressiveRendererTask Declarations
class ProgressiveRendererTask : public Task {
public:
    // ProgressiveRendererTask Public Methods added prc Camera which is used to compute averaged image
    ProgressiveRendererTask(const Scene *sc, Renderer *ren, Camera *c, Camera * prc,
                        ProgressReporter &pr, Sampler *ms, Sample *sam, 
                        bool visIds, int tn, int tc,int sd)
	: reporter(pr)
    {
        scene = sc; renderer = ren; camera = c; progCamera=prc; mainSampler = ms;
        origSample = sam; visualizeObjectIds = visIds; taskNum = tn; taskCount = tc; seed=sd;
    }
    void Run();
private:
    // ProgressiveRendererTask Private Data
    const Scene *scene;
    const Renderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    bool visualizeObjectIds;
    int taskNum, taskCount;
//MC added seed value
	Camera* progCamera;
	int seed;
};



#endif // PBRT_RENDERERS_PROGRESSIVERENDERER_H
