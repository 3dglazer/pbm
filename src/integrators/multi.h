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

// MultiScatteringIntegrator Declarations
class MultiScatteringIntegrator : public VolumeIntegrator {
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
    Spectrum Lmm(const Scene *, const Renderer *, const RayDifferential &ray,
                const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
    Spectrum Lsm(const Scene *, const Renderer *, const RayDifferential &ray,
                 const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;
private:
    // ------  TODO -----
    // Should add light cache structures
    
    // MultiScatteringIntegrator Private Data
    float stepSize;
    int tauSampleOffset, scatterSampleOffset;
};


MultiScatteringIntegrator *CreateMultiScatteringIntegrator(const ParamSet &params);


#endif /* defined(__pbrt__multi__) */