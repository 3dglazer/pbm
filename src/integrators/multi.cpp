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

// Variation of the Code from http://geomalgorithms.com/a07-_distance.html
// References David Eberly, "Distance Methods" in 3D Game Engine Design (2006)
//Seth Teller, line_line_closest_points3d() (2000) cited in the Graphics  Algorithms FAQ (2001)
// Returns closest distance between lines and creates two closest points on these lines.
float ln2lnPP(const Ray &L1,const Ray &L2, float &p1, float &p2)
{
    Vector u = L1.d;
    Vector v = L2.d;
    Vector w = L1.o - L2.o;
    float a = Dot(u,u);         
    float b = Dot(u,v);
    float c = Dot(v,v);         
    float d = Dot(u,w);
    float e = Dot(v,w);
    float D = a*c - b*b;        
    float sc, tc;
    if (D < 0.000001) {          // parallel
        sc = 0.0;
        tc = (b>c ? d/b : e/c);    
    }
    else {
        sc = (b*e - c*d) / D;
        tc = (a*e - b*d) / D;
    }
//    if (p1==0 || p2==0) {
//        Vector   dP = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)
//        return dP.Length();   // return the closest distance
//    }else{
        //p1=new Point(L1.o+sc*u);
        //p2=new Point(L2.o+tc*v);
        //u= Vector(*p1-*p2);
    p1=sc;
    p2=tc;
    return u.Length();
//    }
}
//===================================================================

inline float A(float x, float phi, float h){
    return asinhf((x/h)*sinf(phi));
}

// analytical marginal inverse CDF for inverse squared sampling points vrl, phi is angle between camera vector and vrl vector
inline float invSqrCDF(float r,float phi,float h, float v0, float v1){
    return h*sinhf(Lerp(A(v0, phi, h), A(v1,phi,h), r))/sinf(phi);
}

inline float B(float x,float h){
    return atan(x/h);
}

//generates equiangular samples on the camera ray between u0 and u1 parametric distances, r is random number
inline float invEqAngCDF(float r,float h,float u0,float u1){
    return h*tanf(Lerp(B(u0,h), B(u1,h), r));
}

//both rays have to be normalized + the vrl.mint is beginning of the media and vrl.maxt is the end of interraction. h will contain the smallest distance between the vectors
float sampleVRL(const Ray &r,const Ray &vrl, RNG& rng, float &h){
    float uh; //closest point parameter on camera ray
    float vh; //closest point parameter on vrl
    h=ln2lnPP(r, vrl, uh, vh);
    float phi=acosf(Dot(r.d, vrl.d));
    float v0d= vrl.mint-h;
    float v1d= vrl.maxt-h;
    float rnd=rng.RandomFloat();
    return invSqrCDF(rnd, phi, h, v0d, v1d);
}

//unit vectors expected otherwise wont work
inline float myline2line(const Point &p0,const Vector &v0,const Point &p1,const Vector &v1){
    bool paralel=false;
    Vector tmp=(v0-v1);
    if(tmp.Length()<0.00001) paralel=true;
    if (paralel) {
        tmp=Cross(Vector(p1-p0), v0);
        return tmp.Length();
    }else{
        tmp=Cross(v0,v1);
        return Dot(Vector(p1-p0), tmp)/tmp.Length();
    }
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
    Spectrum Lmm= vrlSamplingBruteForce(t0, t1, scene, renderer, ray, sample, rng, T, arena);
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
    
    Lv=vslSamplingBruteForce(t0, t1, scene, renderer, ray, sample, rng, T, arena);
    //Lv=vslSamplingVIZ(t0, t1, scene, renderer, ray, sample, rng, T, arena);
    
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
                float vrlD=curPath->ray.mint+(curPath->ray.maxt - curPath->ray.mint) * (rng.RandomFloat());
                vrlSample=curPath->ray.o+curPath->ray.d*vrlD;
                Vector wo=Vector(vrlSample - p); //vector to the light
                Ray connectRay= Ray(p, wo, 0.);
                if (scene->IntersectP(connectRay)) {
                    continue; //move on sample is obscured
                }
                vrlTr=curPath->getTransmittance(vrlD);
                Spectrum vrlContrib=curPath->contrib * vrlTr;
                
                float d2=wo.LengthSquared();
                wo=Normalize(wo);
                float pp=vr->p(p, ray.d, wo, ray.time); // phase phunction at current point
                float pvrl=vr->p(vrlSample, curPath->ray.d, -wo, ray.time);
                Spectrum ssVrl= vr->sigma_s(vrlSample, curPath->ray.d, ray.time);
                if (isnan(pp)||isnan(pvrl)) {
                    continue;
                }
                if (!vrlContrib.IsBlack() && !ss.IsBlack() && !ssVrl.IsBlack() && pp!=0 && pvrl!=0) {
                    Lv+=Tr*vrlContrib*ss*pp*pvrl; //vrl ss is counted in the contrib ?? *ssVrl // the square dist term is missing
                    Lv*=1./d2;  //test without phase functions
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
    Spectrum tau;
    Spectrum Lmm=0.f;
    float d;
    //create a VRL
    VolumePath rayTransm=VolumePath(ray);
    rayTransm.ray.mint=t0;
    rayTransm.ray.maxt=t1;
    
    
    //maybe add some constraining criteria like maximum scattering events in the original VRL paper is 16
    for (int evnts=0; evnts < 8; ++evnts) {
        d=renderer->freeFlight(scene, rayTransm.ray, tau, rng); //the tau could be used for multiple scattering
        if (d==-1.){
            break; //scattering did not happened tau should be valid
        }else{
            ray.mint=d; //save event distance for next tracking
            //add the transmittance and distance to the VRL
            rayTransm.dists.push_back(d);
            //phase function doesn't have to be saved can be querried directly from scene->volumeRegion->p(....) returns probability float.
        }
        //possibly terminate tracking if the transmittance is really small
        if (tau.returnOne() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
        }
    }// end of SCATTERING EVENTS
    //the scattering might end before reaching the surface we have to compute the tau
    if (d<t1) {
        tau+=(t1-d)*scene->volumeRegion->sigma_t(ray.o+ray.d*ray.mint, ray.d, ray.time);
    }
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

