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
//#include <cmath>
//#include <math.h>
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
//right now computes single scattering using raymarching to light
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
    Vector wo;
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
Spectrum MultiScatteringIntegrator::Lmm(const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena)  const {
    VolumeRegion *vr = scene->volumeRegion;
    RayDifferential localRay(ray);
    float t0,t1;
    if (!vr || !vr->IntersectP(localRay, &t0, &t1) || (t1-t0) == 0.f) {
        *T = 1.f;
        return 0.f;
    }
#ifdef BRUTE
    Spectrum Lmm= vrlSamplingBruteForce(t0, t1, scene, renderer, ray, sample, rng, T, arena);
#else
    Spectrum Lmm=vrlSamplingPaper1(t0, t1, scene, renderer, ray, sample, rng, T, arena);
#endif
    //Spectrum Lmm= vrlSamplingVIZ(t0, t1, scene, renderer, ray, sample, rng, T, arena);

    return Lmm;
}



//Surface -> Media transport evaluation part, should integrate over VPLs or VSLs
Spectrum MultiScatteringIntegrator::Lsm(const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena)  const{
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);
    VolumeRegion *vr = scene->volumeRegion;
    RayDifferential localRay(ray);
    
    float t0,t1;
    if (!vr || !vr->IntersectP(localRay, &t0, &t1) || (t1-t0) == 0.f) {
        *T = 1.f;
        return 0.f;
    }
#ifdef BRUTE
    Lv=vslSamplingBruteForce(t0, t1, scene, renderer, ray, sample, rng, T, arena);
#else
    //Lv=vslSamplingVIZ(t0, t1, scene, renderer, ray, sample, rng, T, arena);
    Lv=vslSamplingCDF(t0, t1, scene, renderer, ray, sample, rng, T, arena);
#endif
    
    return Lv;
}

MultiScatteringIntegrator *CreateMultiScatteringIntegrator(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
    return new MultiScatteringIntegrator(stepSize);
}

Spectrum MultiScatteringIntegrator::vrlSamplingBruteForce(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const{
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);
    VolumeRegion *vr = scene->volumeRegion;
    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;
    
    // Compute sample patterns for single scattering samples
//    float *lightNum = arena.Alloc<float>(nSamples);
//    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
//    float *lightComp = arena.Alloc<float>(nSamples);
//    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
//    float *lightPos = arena.Alloc<float>(2*nSamples);
//    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
    uint32_t sampOffset = 0;
    
    Point vrlSample;
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
        Spectrum stepTau = vr->tau(tauRay,.5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau);
        
        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }
        
        // Compute single-scattering source term at _p_
        Lv += Tr * vr->Lve(p, w, ray.time); //emmission
        
        Spectrum ss = vr->sigma_s(p, w, ray.time);
        if (!ss.IsBlack() && !this->vpths.empty()) {
            VolumePath *curPath;
            float vrlTr;
            
            for (int vpthIdx=0; vpthIdx<vpths.size(); ++vpthIdx) {
                curPath=vpths[vpthIdx];
                //sample few places on vrl
                for (int did=0; did<10; ++did) {
                    float vrlD=curPath->ray.mint+(curPath->ray.maxt - curPath->ray.mint) * (rng.RandomFloat());
                    vrlSample=curPath->ray.o+curPath->ray.d*vrlD;
                    Vector wo=Normalize(vrlSample - p); //vector to the light
                    Ray connectRay= Ray(p, wo, 0.);
                    if (scene->IntersectP(connectRay)) {
                        continue; //move on sample is obscured
                    }
                    vrlTr=curPath->getTransmittance(vrlD);
                    Spectrum vrlContrib=curPath->contrib * vrlTr;
                    float d2=DistanceSquared(vrlSample, p);
                    float pp=vr->p(p, w, wo, ray.time); // phase phunction at current point
                    float pvrl=vr->p(vrlSample, curPath->ray.d, -wo, ray.time);
                    //Spectrum ssVrl= vr->sigma_s(vrlSample, curPath->ray.d, ray.time);
                    if (isnan(pp)||isnan(pvrl)) {
                        continue;
                    }
                    //&& !ssVrl.IsBlack()&& !ssVrl.IsBlack()
                    if (!vrlContrib.IsBlack() && !ss.IsBlack() && pp!=0 && pvrl!=0) {
                        Lv+=Tr*vrlContrib*pp*pvrl*renderer->Transmittance(scene, RayDifferential(connectRay), NULL, rng, arena)*1./d2*ss;//*ssVrl; //vrl ss is counted in the contrib ?? *ssVrl // the square dist term is missing
                    }
                }
            }
             
//            int nLights = scene->lights.size();
//            int ln = min(Floor2Int(lightNum[sampOffset] * nLights),
//                         nLights-1);
//            Light *light = scene->lights[ln];
            // Add contribution of _light_ due to scattering at _p_
//            float pdf;
//            VisibilityTester vis;
//            Vector wo;
//            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
//                           lightPos[2*sampOffset+1]);
//            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            
//            if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {
//                Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);
//                Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * Ld * float(nLights) /
//                pdf;
//            }
        }
        ++sampOffset;
    }
    *T = Tr;
    return Lv * step;
}

Spectrum MultiScatteringIntegrator::vrlSamplingVIZ(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const{
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);
    VolumeRegion *vr = scene->volumeRegion;
    // Prepare for volume integration stepping
    Point vrlSample;
    
    float vrlRadius=0.0001;
    if ( !this->vpths.empty()) {
            VolumePath *curPath;
            for (int vpthIdx=0; vpthIdx<vpths.size(); ++vpthIdx) {
                curPath=vpths[vpthIdx];
                float p0,p1;
                float dist=myline2line(ray.o, Normalize(ray.d), curPath->ray.o, Normalize(curPath->ray.d));
                //float dist=ln2lnPP(ray,curPath->ray,p0,p1);
                if ((dist<0?-dist:dist)<vrlRadius) {
                    Lv+=0.001;
                }
            }

        }
    return Lv ;
}

Spectrum MultiScatteringIntegrator::vslSamplingVIZ(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const{
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);
    VolumeRegion *vr = scene->volumeRegion;
    // Prepare for volume integration stepping
    Point vrlSample;
    
    float vrlRadius=1.3;
    if ( !this->vsls.empty()) {
        VirtualSphericalLight *curVsl;
        for (int vpthIdx=0; vpthIdx<vsls.size(); ++vpthIdx) {
            curVsl=vsls[vpthIdx];
            float p0,p1;
            
            float dist=(Cross( Vector(ray.o-curVsl->p),ray.d)).Length()/ray.d.Length();
            if (dist<curVsl->radius) {
                Lv+=1.3;
            }
        }
        
    }
    return Lv ;
}

Spectrum MultiScatteringIntegrator::vrlSamplingPaper1(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const{
    
#ifndef OLDIMPL
    
    VolumeRegion *vr = scene->volumeRegion;
    Vector w = -ray.d;
    RayDifferential lray(ray);
    //sample every vrl nSamples times
    const int nSamples=100;
    const float stepSize=1./(float)nSamples;
    float pdfVals[nSamples];
    float rayPoint=0;
    const int camSamples=16;
    //first save transmitance to VolumePath take from particleShooter
    //maybe add some constraining criteria like maximum scattering events in the original VRL paper is 16
    Spectrum tau(0.);
    Spectrum maxtau;
    float d;
    float dmax=-INFINITY;
    VolumePath rayTransmCache(ray);
    //rayTransmCache.ray.d=-ray.d;
    rayTransmCache.ray.mint=t0;
    rayTransmCache.ray.maxt=t1;
    Spectrum vslsContrib;
    Spectrum L(0.);
    Spectrum alpha(1.);
    for (int evnts=0; evnts < camSamples; ++evnts) {
        d=renderer->freeFlight(scene, lray, maxtau, rng); //the tau could be used for multiple scattering
        if (d==-1.){
            rayTransmCache.dists.push_back(ray.maxt);
            continue; //scattering did not happened tau should be valid
        }else{
            if (d>dmax) {
                dmax=d;
                tau=maxtau;
            }
            rayTransmCache.dists.push_back(d);
        }
    }// end of SCATTERING EVENTS
    
    //if full transmitance needed
    if (dmax<t1) {
        alpha*=renderer->Transmittance(scene, ray, NULL, rng, arena); //attnuation to the vsl
        //tau+=(t1-d)*scene->volumeRegion->sigma_t(ray.o+ray.d*t1, ray.d, ray.time);
    }
    
    *T = alpha;
    
    
    if (vsls.empty()) {
        return L;
    }
    
    
    //for now fixed point in the middle of the vrl is used 
    float vrlSampleDist;
    Ray mray=rayTransmCache.ray;
    
    //nead point to line
    Point minP=mray.o+mray.d*mray.mint;
    Point maxP=mray.o+mray.d*mray.maxt;
    Point p;
    Point vpthPoint;
    //iterating over Volume paths sigle sample on every vpl
    for (uint32_t i = 0; i < this->vpths.size(); ++i) {
        //printf("\n======iterating over vsls========");
        VolumePath *currVrl = vpths[i];
        //for (int ksicht=0; ksicht<10; ++ksicht) {
        vrlSampleDist=currVrl->ray.mint+(currVrl->ray.maxt-currVrl->ray.mint)*rng.RandomFloat();
        //vrlSampleDist=currVrl->ray.mint;
        Point vrlPoint=currVrl->ray.o+currVrl->ray.d*(vrlSampleDist);
        Spectrum ssVrl= vr->sigma_s(vrlPoint, currVrl->ray.d, ray.time);
        float sParam=0.; //parametric distance to vector beginning
        float h=p2Ray(vrlPoint, minP,maxP,sParam);
        //printf("\nh=%f; sParam=%f;\n",h,sParam);
        float v0=mray.mint-sParam;
        float v1=mray.maxt-sParam;
        int smpl=0;
        // generate analytical marginal pdf for vrl Sampling using 10 equiAngular samples on VRL
        for (float s=0; s<=1.; s+=stepSize) {
            //generate 10 eq angular samples
            rayPoint=invEqAngCDF(s, h, v0, v1);
            rayPoint+=sParam; //move the sample to ray origin

            p=minP+mray.d*rayPoint;
            Vector wi = Normalize(vrlPoint - p);
            float pp=vr->p(p, w, -wi, ray.time); // phase phunction at current pointat current point
            float ppVrl=vr->p(vrlPoint, currVrl->ray.d, wi, ray.time);
            float G = pp * ppVrl;
            pdfVals[smpl]=G;
            smpl++;
        }
        
        //sample according to the pdf distribution
        float ret,pdf;
        Distribution1D distrib=Distribution1D(&pdfVals[0], nSamples);
        // sample ten points on vrl acording to the inverse squared cdf
        for (int vrlS=0; vrlS<10; ++vrlS) {
            ret=distrib.SampleContinuous(rng.RandomFloat(),&pdf); //what is the pdf and what is the ret
            //generate 10 eq angular samples
            rayPoint=invEqAngCDF(ret, h, v0, v1);
            rayPoint+=sParam; //move the sample to ray origin
            
            p=minP+mray.d*rayPoint;
            float d2 = DistanceSquared(p, vrlPoint);
            Vector wi = Normalize(vrlPoint - p);
            RayDifferential connectRay(p, wi, ray, NULL, sqrtf(d2));
            if (scene->IntersectP(connectRay)) { //move on the next light if the light is obscured
                continue;
            }
            float pp=vr->p(p, w, -wi, ray.time); // phase phunction at current pointat current point
            float ppVrl=vr->p(vrlPoint, currVrl->ray.d, wi, ray.time);
            float G = pp * ppVrl/d2;
            G = (G<10000.)?G:10000.;
            if (G == 0.f) continue;
            Spectrum Llight = G * currVrl->contrib*currVrl->getTransmittance(vrlSampleDist);
            //Spectrum Llight = G * vl->pathContrib; //neuvazuju brdf jen pro zkousku
            Spectrum ss = vr->sigma_s(p, w, ray.time);
            
            if (ss.IsBlack()||ssVrl.IsBlack()) {
                continue;
            }
            Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, arena);
            float rdstnc=rayPoint+mray.mint;
            float Tr =rayTransmCache.getTransmittance(rdstnc);
            //possibly compute single scattering term
            L += Tr * vr->Lve(p, w, ray.time); //emmission
            if (Tr==0) {
                if (!alpha.IsBlack()) {
                    Tr=Lerp(rdstnc, 1., alpha.y());
                }
                continue;
            }
            if (!ss.IsBlack())
                L += Tr*ss*ssVrl*Llight;
           // }
        }
    }

    return L;
    
    
    
#else
    
    
    
    
    
    
    
    
    Spectrum tau;
    Spectrum Lmm=0.f;
    float d;
    //create a VRL
    VolumePath rayTransm=VolumePath(ray);
    rayTransm.ray.mint=t0;
    rayTransm.ray.maxt=t1;
    
    
    Spectrum maxtau;
    float dmax=-INFINITY;
    for (int evnts=0; evnts < 16; ++evnts) {
        d=renderer->freeFlight(scene, ray, maxtau, rng); //the tau could be used for multiple scattering
        if (d==-1.){
            rayTransm.dists.push_back(ray.maxt);
            continue; //scattering did not happened tau should be valid
        }else{
            if (d>dmax) {
                dmax=d;
                tau=maxtau;
            }
            //ray.mint=d; //save event distance for next tracking
            //add the transmittance and distance to the VRL
            rayTransm.dists.push_back(d);
            //phase function doesn't have to be saved can be querried directly from scene->volumeRegion->p(....) returns probability float.
        }
    }// end of SCATTERING EVENTS
    //if the tracking ended before reaching the end point of media we have to add the optical thickness of the rest of the media to it, to correctly attenuate the vsls.
    if (dmax<t1) {
        *T=renderer->Transmittance(scene, ray, NULL, rng, arena); //attnuation to the vsl
        //tau+=(t1-d)*scene->volumeRegion->sigma_t(ray.o+ray.d*t1, ray.d, ray.time);
    }
    std::sort((rayTransm.dists.begin()), (rayTransm.dists.end()));
    //computed transmittance to be returned
    *T=Exp(-tau);
    
    VolumePath *curPath;
    //iterating over Volume paths
    for (int vpthIdx=0; vpthIdx<vpths.size(); ++vpthIdx) {
        curPath=vpths[vpthIdx];
        //sample point on this vrl
        float h; // will contain the smallest destance between VRL and ray
        float vrlDist=sampleVRL(ray, curPath->ray, rng,h);
        // for now 10 samples will be generated along camera ray
        float rayParam;
        float vrlTr, rayTr, vrlP, rayP;
        Point vrlSample=curPath->ray.o + curPath->ray.d * vrlDist;
        Point raySample;
        Vector sampleVect;
        // TODO add scattering terms!!!!!!!
        for (int camRaySampleId=0; camRaySampleId<10; ++camRaySampleId) {
            //suboptimal sampling in heterogenous media!!!
            rayParam=invEqAngCDF(rng.RandomFloat(),h, ray.mint, ray.maxt);
            vrlTr=curPath->getTransmittance(vrlDist); //gets VRL transmittance
            rayTr=rayTransm.getTransmittance(rayParam); // gets transmittance along ray
            raySample=ray.o + ray.d * rayParam;
            sampleVect= Vector(vrlSample-raySample);
            rayP= scene->volumeRegion->p(raySample, ray.d, sampleVect, ray.time);
            vrlP= scene->volumeRegion->p(vrlSample, curPath->ray.d, -sampleVect, ray.time);
            if (vrlTr==0||rayTr==0||rayP==0||vrlP==0||isnan(rayP)||isnan(vrlP)) {
                continue;
            }
            Lmm+=(curPath->contrib * vrlTr * rayTr) * (rayP * vrlP);
        }
    }
    return Lmm;
#endif
    
    
}

Spectrum MultiScatteringIntegrator::vslSamplingBruteForce(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const{
    VolumeRegion *vr = scene->volumeRegion;
    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;
    
    uint32_t sampOffset = 0;
    
    Point vrlSample;    
    Spectrum Lv=0.f;
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
        Spectrum stepTau = vr->tau(tauRay,.5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau);
        
        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }
        
        // Compute single-scattering source term at _p_
        Lv += Tr * vr->Lve(p, w, ray.time); //emmission
        
        Spectrum ss = vr->sigma_s(p, w, ray.time);
        if (!ss.IsBlack() && !this->vsls.empty()) {
            //iterating over Volume paths sigle sample on every vpl
            for (uint32_t i = 0; i < this->vsls.size(); ++i) {
                //printf("\n======iterating over vsls========");
                VirtualSphericalLight *vl = vsls[i];
                // Compute virtual light's tentative contribution _Llight_
                float d2 = DistanceSquared(p, vl->p);
                Vector wi = Normalize(vl->p - p);
                RayDifferential connectRay(p, wi, ray, NULL, sqrtf(d2) * (1.f - vl->rayEpsilon));
                if (scene->IntersectP(connectRay)) { //move on the next light if the light is obscured
                    continue;
                }
                float pp=vr->p(p, w, -wi, ray.time); // phase phunction at current point
                float G = pp * AbsDot(wi, vl->n) / d2;
                G = (G<10000.)?G:10000.;
                Spectrum f = vl->bsdf->f(-wi, vl->i); // is -wi correct??
//                //Spectrum f = bsdf->f(wo, wi);
                if (G == 0.f || f.IsBlack()) continue;
                Spectrum Llight = f * G * vl->pathContrib;
                //Spectrum Llight = G * vl->pathContrib; //neuvazuju brdf jen pro zkousku
                
                Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, arena);
                if (!ss.IsBlack())
                    Lv += Tr*ss*Llight;
            }
        }
        ++sampOffset;
    }
    *T = Tr;
    return Lv * step;
}

Spectrum MultiScatteringIntegrator::vslSamplingCDF(float t0,float t1,const Scene *scene, const ProgressiveRenderer * renderer, const RayDifferential &ray,const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const{
    VolumeRegion *vr = scene->volumeRegion;
    Vector w = -ray.d;
    RayDifferential lray(ray);
    //sample every vrl nSamples times
    const int nSamples=10;
    const float stepSize=1./(float)nSamples;
    float pdfVals[nSamples];
    float rayPoint=0;
    const int camSamples=10;
    //first save transmitance to VolumePath take from particleShooter
    //maybe add some constraining criteria like maximum scattering events in the original VRL paper is 16
    Spectrum tau(0.);
    Spectrum maxtau;
    float d;
    float dmax=-INFINITY;
    VolumePath rayTransmCache(ray);
    //rayTransmCache.ray.d=-ray.d;
    rayTransmCache.ray.mint=t0;
    rayTransmCache.ray.maxt=t1;
    Spectrum vslsContrib;
    Spectrum L(0.);
    Spectrum alpha(1.);
    for (int evnts=0; evnts < camSamples; ++evnts) {
        d=renderer->freeFlight(scene, lray, maxtau, rng); //the tau could be used for multiple scattering
        if (d==-1.){
            rayTransmCache.dists.push_back(lray.maxt);
            continue; //scattering did not happened tau should be valid
        }else{
            if (d>dmax) {
                dmax=d;
                tau=maxtau;
            }
            rayTransmCache.dists.push_back(d);
        }
    }// end of SCATTERING EVENTS
    
    //if full transmitance needed
    if (dmax<t1) {
        alpha*=renderer->Transmittance(scene, ray, NULL, rng, arena); //attnuation to the vsl
        //tau+=(t1-d)*scene->volumeRegion->sigma_t(ray.o+ray.d*t1, ray.d, ray.time);
    }
    
    *T = alpha;
    
    
    if (vsls.empty()) {
        return L;
    }
    
    
    Ray mray=rayTransmCache.ray;
    
    //nead point to line
    Point minP=mray.o+mray.d*mray.mint;
    Point maxP=mray.o+mray.d*mray.maxt;
    Point p;
    Point vpthPoint;
    //iterating over Volume paths sigle sample on every vpl
    for (uint32_t i = 0; i < this->vsls.size(); ++i) {
        //printf("\n======iterating over vsls========");
        VirtualSphericalLight *vl = vsls[i];
        
        
        float sParam=0.; //parametric distance to vector beginning
        float h=p2Ray(vl->p, minP,maxP,sParam);
        //printf("\nh=%f; sParam=%f;\n",h,sParam);
        float v0=mray.mint-sParam;
        float v1=mray.maxt-sParam;
        int smpl=0;
        // generate analytical marginal pdf for vrl Sampling using 10 equiAngular samples on VRL
        for (float s=0; s<=1.; s+=stepSize) {
            //generate 10 eq angular samples
            rayPoint=invEqAngCDF(s, h, v0, v1);
            rayPoint+=sParam; //move the sample to ray origin

            p=minP+mray.d*rayPoint;
            Vector wi = Normalize(vl->p - p);
            float pp=vr->p(p, w, -wi, ray.time); // phase phunction at current pointat current point
            float G = pp * AbsDot(wi, vl->n);
            //G = (G<10000.)?G:10000.;
            Spectrum f = vl->bsdf->f(-wi, vl->i);
            if (G == 0.f || f.IsBlack()){
                pdfVals[smpl]=0;
            }else{
                pdfVals[smpl]=G*f.y();
            }
            smpl++;
        }
        
        //sample according to the pdf distribution
        float ret,pdf;
        Distribution1D distrib=Distribution1D(&pdfVals[0], nSamples);
        // sample ten points on vrl acording to the inverse squared cdf
        for (int vrlS=0; vrlS<10; ++vrlS) {
            ret=distrib.SampleContinuous(rng.RandomFloat(),&pdf); //what is the pdf and what is the ret
            //generate 10 eq angular samples
            rayPoint=invEqAngCDF(ret, h, v0, v1);
            rayPoint+=sParam; //move the sample to ray origin
            
            p=minP+mray.d*rayPoint;
            float d2 = DistanceSquared(p, vl->p);
            Vector wi = Normalize(vl->p - p);
            RayDifferential connectRay(p, wi, ray, NULL, sqrtf(d2) * (1.f - vl->rayEpsilon));
            if (scene->IntersectP(connectRay)) { //move on the next light if the light is obscured
                continue;
            }
            float pp=vr->p(p, w, -wi, ray.time); // phase phunction at current point
            float G = pp * AbsDot(wi, vl->n) / d2;
            G = (G<10000.)?G:10000.;
            Spectrum f = vl->bsdf->f(-wi, vl->i); // is -wi correct??
            //                //Spectrum f = bsdf->f(wo, wi);
            if (G == 0.f || f.IsBlack()) continue;
             Spectrum Llight = f * G * vl->pathContrib;
            //Spectrum Llight = G * vl->pathContrib; //neuvazuju brdf jen pro zkousku
            Spectrum ss = vr->sigma_s(p, w, ray.time);
            Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, arena);
            float rdstnc=rayPoint+mray.mint;
            float Tr =rayTransmCache.getTransmittance(rdstnc);
            if (Tr==0) {
                if (!alpha.IsBlack()) {
                    Tr=Lerp(rdstnc, 1., alpha.y());
                }
                continue;
            }
            if (!ss.IsBlack())
                L += Tr*ss*Llight;
            
        }
    }
//
//    Spectrum Tr(0.);
//    const Spectrum one(1.);
//    alpha*=Exp(-tau); //transmittance from optical thickness
//    for (float raydst=0.; raydst<=1.0; raydst+=0.03) {
//        float raylength=(rayTransmCache.ray.maxt-rayTransmCache.ray.mint);
//        Tr =one*rayTransmCache.getTransmittance(raydst*raylength);
//        Point p=rayTransmCache.ray.o+rayTransmCache.ray.d*raylength*raydst;
//        
//        // Compute single-scattering source term at _p_
//        //Lv += Tr * vr->Lve(p, w, ray.time); //emmission
//        Spectrum ss = vr->sigma_s(p, w, ray.time);
//        if (!ss.IsBlack() && !this->vsls.empty()) {
//            //iterating over Volume paths sigle sample on every vpl
//            for (uint32_t i = 0; i < this->vsls.size(); ++i) {
//                //printf("\n======iterating over vsls========");
//                VirtualSphericalLight *vl = vsls[i];
//                // Compute virtual light's tentative contribution _Llight_
//                float d2 = DistanceSquared(p, vl->p);
//                Vector wi = Normalize(vl->p - p);
//                RayDifferential connectRay(p, wi, ray, NULL, sqrtf(d2) * (1.f - vl->rayEpsilon));
//                if (scene->IntersectP(connectRay)) { //move on the next light if the light is obscured
//                    continue;
//                }
//                float pp=vr->p(p, w, -wi, ray.time); // phase phunction at current point
//                float G = pp * AbsDot(wi, vl->n) / d2;
//                G = (G<10000.)?G:10000.;
//                Spectrum f = vl->bsdf->f(-wi, vl->i); // is -wi correct??
//                //                //Spectrum f = bsdf->f(wo, wi);
//                if (G == 0.f || f.IsBlack()) continue;
//                Spectrum Llight = f * G * vl->pathContrib;
//                //Spectrum Llight = G * vl->pathContrib; //neuvazuju brdf jen pro zkousku
//                
//                Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, arena);
//                if (!ss.IsBlack())
//                    L += Tr*ss*Llight;
//            }
//        }
//    }
    
    return L;
}


