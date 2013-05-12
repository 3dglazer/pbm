//
//  multi.h
//  pbrt
//
//  Created by Zdenek Glazer on 03.05.13.
//
//

#ifndef __pbrt__multi__
#define __pbrt__multi__

#include <iostream>
// integrators/single.h*
#include "volume.h"
#include "integrator.h"
#include <math.h>
#include <stdlib.h>
#include "samplingmethods.h"
// MultiScatteringIntegrator Declarations
class MultiScatteringIntegrator : public ProgressiveVolumeIntegrator {
public:
    // MultiScatteringIntegrator Public Methods
    MultiScatteringIntegrator(float ss) { stepSize = ss; }
    Spectrum Transmittance(const Scene *, const Renderer *,
                           const RayDifferential &ray, const Sample *sample, RNG &rng,
                           MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample,
                        const Scene *scene);
    Spectrum Li(const Scene *, const Renderer *, const RayDifferential &ray,
                const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    float freeFlight(const Scene *scene, const Ray &r,Spectrum& tau,const RNG &rng) const;
    Spectrum Lmm(const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena)  const;
    Spectrum Lsm(const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena)  const;
private:
    // ------  TODO -----
    // Should add light cache structures
    
    // MultiScatteringIntegrator Private Data
    // VRL sampling methods
    Spectrum vrlSamplingPaper1(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    Spectrum vrlSamplingBruteForce(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    Spectrum vrlSamplingVIZ(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    // VPL sampling methods
    Spectrum vslSamplingBruteForce(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    Spectrum vslSamplingVIZ(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    Spectrum vslSamplingCDF(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    float stepSize;
    int tauSampleOffset, scatterSampleOffset;
};


MultiScatteringIntegrator *CreateMultiScatteringIntegrator(const ParamSet &params);


#endif /* defined(__pbrt__multi__) */