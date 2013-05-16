
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

 */


// integrators/igi.cpp*
#include "stdafx.h"
#include "integrators/igi.h"
#define BRUTE
// IGIIntegrator Method Definitions
IGIIntegrator::~IGIIntegrator() {
	
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
	
	if (dump==false) {
		return;
	}
	igiLocalArena.FreeAll();
	printf("file in igi is: %s",filename.c_str());
}


void IGIIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                   const Scene *scene) {
    // Allocate and request samples for sampling all lights
    uint32_t nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }
    vlSetOffset = sample->Add1D(1);

    if (sampler) nGatherSamples = sampler->RoundSize(nGatherSamples);
    gatherSampleOffset = BSDFSampleOffsets(nGatherSamples, sample);
}

void IGIIntegrator::Preprocess(const Scene *scene, const Camera *camera,
                               const Renderer *renderer) {
	
}

Spectrum IGIIntegrator::Lms(const Scene *scene, const ProgressiveRenderer *renderer,
                            const RayDifferential &ray, const Intersection &isect,
                            const Sample *sample, RNG &rng, MemoryArena &localArena) const {
    //should be called in every raycast in lss for every intersection point
#ifdef BRUTE
    Spectrum ms=this->sampleVRLBruteForce(scene,renderer,ray,isect,sample,rng,localArena);
#else
   // Spectrum ms=this->Li(scene, renderer, ray, isect, sample, rng, localArena);
    Spectrum ms=this->sampleVRLCDF(scene,renderer,ray,isect,sample,rng,localArena);
#endif
   return ms;
    return NULL;
}

Spectrum IGIIntegrator::Lss(const Scene *scene, const ProgressiveRenderer *renderer,
                            const RayDifferential &ray, const Intersection &isect,
                            const Sample *sample, RNG &rng, MemoryArena &localArena) const {
    //this code should be same as in the IGI integrator
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);
    
    Spectrum ss;
    return ss;
}


Spectrum IGIIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &localArena) const {
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);
    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, localArena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    //sample all lights uniformly for dirrect lighting
    L += UniformSampleAllLights(scene, renderer, localArena, p, n,
                    wo, isect.rayEpsilon, ray.time, bsdf, sample, rng,
                    lightSampleOffsets, bsdfSampleOffsets);
    for (uint32_t i = 0; i < this->vsls.size(); ++i) {
        if (rng.RandomFloat()>0.5) {
            continue;
        }
		//printf("\n======iterating over vsls========");
        VirtualSphericalLight *vl = vsls[i];
        // Compute virtual light's tentative contribution _Llight_
        float d2 = DistanceSquared(p, vl->p);
        Vector wi = Normalize(vl->p - p);
        float G = AbsDot(wi, n) * AbsDot(wi, vl->n) / d2;
        //G = min(G, gLimit);
        Spectrum f = bsdf->f(wo, wi);
        if (G == 0.f || f.IsBlack()) continue;
        Spectrum Llight = f * G * vl->pathContrib;
        
        // Possibly skip virtual light shadow ray with Russian roulette
        if (Llight.y() < rrThreshold) {
            float continueProbability = .1f;
            if (rng.RandomFloat() > continueProbability)
                continue;
            Llight /= continueProbability;
        }
        
        RayDifferential connectRay(p, wi, ray, isect.rayEpsilon,
                                   sqrtf(d2) * (1.f - vl->rayEpsilon));
        if (scene->IntersectP(connectRay)) continue;
        Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, localArena);

        // Add contribution from _VirtualLight_ _vl_
        L += Llight;
    }
    if (ray.depth + 1 < maxSpecularDepth) {
        Vector wi;
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             localArena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              localArena);
    }
    return L;
}



Spectrum IGIIntegrator::sampleVRLCDF(const Scene *scene, const Renderer *renderer, const RayDifferential &ray, const Intersection &isect, const Sample *sample, RNG &rng, MemoryArena &localArena) const {
    Assert(!ray.HasNaNs());
    VolumeRegion* vr=scene->volumeRegion;
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);
    
    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, localArena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    VolumePath* curVpth;
    Point vpthPoint;
    Spectrum vpthContrib;
    
    //sample every vrl 5 times
    const int nSamples=10;
    const float stepSize=1./(float)nSamples;
    float pdfVals[nSamples];
    float rayPoint=0;
    //evaluate every vrl
    for (int i=0; i<vpths.size(); ++i) {
        curVpth=vpths[i];

        //nead point to line
        Point minP=curVpth->ray.o+curVpth->ray.d*curVpth->ray.mint;
        Point maxP=curVpth->ray.o+curVpth->ray.d*curVpth->ray.maxt;
        float sParam=0.; //parametric distance to vector beginning
        Ray r;
        float h=p2Ray(p, minP,maxP,sParam);
        //printf("\nh=%f; sParam=%f;\n",h,sParam);
        float v0=curVpth->ray.mint-sParam;
        float v1=curVpth->ray.maxt-sParam;
        int smpl=0;
        // generate analytical marginal pdf for vrl Sampling using 10 equiAngular samples on VRL
        for (float s=0; s<=1.; s+=stepSize) {
            //generate 10 eq angular samples
            rayPoint=invEqAngCDF(s, h, v0, v1);
            rayPoint+=sParam; //move the sample to ray origin
            //printf("rayPoint=%f; iter=%f; ray.mint=%f; h=%f; sParam=%f;\n",rayPoint,s,curVpth->ray.mint,h,sParam);
            vpthPoint=minP+curVpth->ray.d*rayPoint;
            Vector wi=Normalize(vpthPoint-p);
            //float d2 =DistanceSquared(vpthPoint, p);
            //RayDifferential connectRay(p,wi,NULL, isect.rayEpsilon,NULL);
            RayDifferential connectRay;
            connectRay.o=p;
            connectRay.d=wi;
            connectRay.mint=isect.rayEpsilon;
            float pp=vr->p(vpthPoint, curVpth->ray.d, -wi, ray.time); // phase phunction at current point
            float G = pp * AbsDot(wi, n);// / d2;disnace is already accounted with in the sampling scheme
            Spectrum f = bsdf->f(wo, wi);
            if (G == 0.f || f.IsBlack()){
                pdfVals[smpl]=0;
            }else{
                pdfVals[smpl]=G*f.y();
            }
            smpl++;
        }
        float ret,pdf;
        Distribution1D distrib=Distribution1D(&pdfVals[0], nSamples);
        // sample ten points on vrl acording to the inverse squared cdf
        for (int vrlS=0; vrlS<10; ++vrlS) {
            ret=distrib.SampleContinuous(rng.RandomFloat(),&pdf); //what is the pdf and what is the ret
            rayPoint=invEqAngCDF(ret, h, v0, v1);
            rayPoint+=sParam; //move the sample to ray origin
            vpthPoint=minP+curVpth->ray.d*rayPoint;
            //printf("\ncdf Sampled rayPoint=%f; ret=%f;\n",rayPoint,ret);
            Vector wi=Normalize(vpthPoint-p);
            float d2 =DistanceSquared(vpthPoint, p);
            RayDifferential connectRay;
            connectRay.o=p;
            connectRay.d=wi;
            connectRay.mint=isect.rayEpsilon;
            //if it's ocluded continue
            if (scene->IntersectP(connectRay)) {
                continue;
            }
            float pp=vr->p(vpthPoint, curVpth->ray.d, -wi, ray.time); // phase phunction at current point
            float G = pp * AbsDot(wi, n) / d2;
            //G = (G<10.)?G:10.;
            G = (G<100.)?G:100.;
            Spectrum f = bsdf->f(wo, wi);
            if (G == 0.f || f.IsBlack()) continue;
            //weight contribution with vrl transmittance and transmittance between surface point and sample point on vrl
            vpthContrib=curVpth->contrib*curVpth->getTransmittance(rayPoint)*renderer->Transmittance(scene, connectRay, NULL, rng, localArena);
            //weight the contribution with scattering coeficient of the media in the sample point the vector here is not needed
            vpthContrib*=vr->sigma_s(vpthPoint, wi, ray.time);
            //weight the contribution with cos(theta) plus inverse squared distance pluss phase function in the vrl sample point
            vpthContrib*=G;
            //weight the contribution with surface brdf
            L+=vpthContrib*f;
        }
    }
    return L;
}

Spectrum IGIIntegrator::sampleVRLBruteForce(const Scene *scene, const Renderer *renderer,
                                              const RayDifferential &ray, const Intersection &isect,
                                              const Sample *sample, RNG &rng, MemoryArena &localArena) const {
    Assert(!ray.HasNaNs());
    VolumeRegion* vr=scene->volumeRegion;
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);
    
    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, localArena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    VolumePath* curVpth;
    Point vpthPoint;
    Spectrum vpthContrib;
    //evaluate every vrl
    for (int i=0; i<vpths.size(); ++i) {
        if (rng.RandomFloat()>0.5) {
            continue;
        }
        curVpth=vpths[i];
        //sample every vrl 5 times
        for (int s=0; s<10; ++s) {
            float d=rng.RandomFloat();
            float rayPoint=curVpth->ray.mint+(curVpth->ray.maxt-curVpth->ray.mint)*d;
            vpthPoint=curVpth->ray.o+curVpth->ray.d*rayPoint;
            Vector wi=Normalize(vpthPoint-p);
            float d2 =DistanceSquared(vpthPoint, p);
            RayDifferential connectRay;
            connectRay.o=p;
            connectRay.d=wi;
            connectRay.mint=isect.rayEpsilon;
            //if it's ocluded continue
            if (scene->IntersectP(connectRay)) {
                continue;
            }
            
            float pp=vr->p(vpthPoint, curVpth->ray.d, -wi, ray.time); // phase phunction at current point
            float G = pp * AbsDot(wi, n) / d2;
            //float G= 1.f/d2;
            G = (G<10.)?G:10.;
            Spectrum f = bsdf->f(wo, wi);
            if (G == 0.f || f.IsBlack()) continue;
            //weight contribution with vrl transmittance and transmittance between surface point and sample point on vrl
            vpthContrib=curVpth->contrib*curVpth->getTransmittance(rayPoint)*renderer->Transmittance(scene, connectRay, NULL, rng, localArena);
            //weight the contribution with scattering coeficient of the media in the sample point the vector here is not needed 
            vpthContrib*=vr->sigma_s(vpthPoint, wi, ray.time);
            //weight the contribution with cos(theta) plus inverse squared distance pluss phase function in the vrl sample point
            vpthContrib*=G;
            //weight the contribution with surface brdf
            L+=vpthContrib*f;

        }
    }
    return L;
}

IGIIntegrator *CreateIGISurfaceIntegrator(const ParamSet &params) {
    int nLightPaths = params.FindOneInt("nlights", 64);
    if (PbrtOptions.quickRender) nLightPaths = max(1, nLightPaths / 4);
    int nLightSets = params.FindOneInt("nsets", 4);
    float rrThresh = params.FindOneFloat("rrthreshold", .0001f);
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float glimit = params.FindOneFloat("glimit", 10.f);
    int gatherSamples = params.FindOneInt("gathersamples", 16);
	//MC stores image filename from params
	string filename = params.FindOneString("filename", PbrtOptions.imageFile);
	bool dump=params.FindOneBool("dump", false);
    return new IGIIntegrator(nLightPaths, nLightSets, rrThresh,
                             maxDepth, glimit, gatherSamples,dump,filename);
	
	
}


