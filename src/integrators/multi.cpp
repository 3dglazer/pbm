//
//  multi.cpp
//  pbrt
//
//  Created by Zdenek Glazer on 03.05.13.
//
//
#include "stdafx.h"
#include "integrators/multi.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"

// MultiScatteringIntegrator Method Definitions
void MultiScatteringIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                                const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

// Integrates volume thickness and returns transmitance
Spectrum MultiScatteringIntegrator::Transmittance(const Scene *scene,
                                                   const Renderer *renderer, const RayDifferential &ray,
                                                   const Sample *sample, RNG &rng, MemoryArena &arena) const {
    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}


//Compatibility solution for pbrt should evaluate both Lmm and Lsm + transmittance to the closest visible surface.
Spectrum MultiScatteringIntegrator::Li(const Scene *scene, const Renderer *renderer,
                                        const RayDifferential &ray, const Sample *sample, RNG &rng,
                                        Spectrum *T, MemoryArena &arena) const {
    VolumeRegion *vr = scene->volumeRegion;
    float t0, t1;
    if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
        *T = 1.f;
        return 0.f;
    }
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);
    
    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;
    
    // Compute sample patterns for single scattering samples
    float *lightNum = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
    float *lightComp = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
    float *lightPos = arena.Alloc<float>(2*nSamples);
    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
    uint32_t sampOffset = 0;
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
        Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau);
        
        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }
        
        // Compute single-scattering source term at _p_
        Lv += Tr * vr->Lve(p, w, ray.time);
        Spectrum ss = vr->sigma_s(p, w, ray.time);
        if (!ss.IsBlack() && scene->lights.size() > 0) {
            int nLights = scene->lights.size();
            int ln = min(Floor2Int(lightNum[sampOffset] * nLights),
                         nLights-1);
            Light *light = scene->lights[ln];
            // Add contribution of _light_ due to scattering at _p_
            float pdf;
            VisibilityTester vis;
            Vector wo;
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            
            if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {
                Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);
                Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * Ld * float(nLights) /
                pdf;
            }
        }
        ++sampOffset;
    }
    *T = Tr;
    return Lv * step;
}

//MC added volume tracking
float MultiScatteringIntegrator::freeFlight(const Scene *scene, const Ray &r,Spectrum& tau,const RNG &rng) const {
    if (!scene->volumeRegion) return -1.0;
    return scene->volumeRegion->freeFlight(r, tau, rng);
}


//Media -> Media integration part integration over VRLs
Spectrum MultiScatteringIntegrator::Lmm(const Scene *, const ProgressiveRenderer *, const RayDifferential &ray,
                                        const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena)  const {
    Spectrum Tr(1.f);
    return Tr;
}

//Surface -> Media transport evaluation part, should integrate over VPLs or VSLs
Spectrum MultiScatteringIntegrator::Lsm(const Scene *, const ProgressiveRenderer *, const RayDifferential &ray,
                                        const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena)  const{
    Spectrum Tr(1.f);
    return Tr;
}
MultiScatteringIntegrator *CreateMultiScatteringIntegrator(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
    return new MultiScatteringIntegrator(stepSize);
}
