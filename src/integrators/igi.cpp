
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
#include "scene.h"
#include "montecarlo.h"
#include "progressreporter.h"
#include "sampler.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"


//MC tests with volumetric photon mapping
#include "photonDisc.h"
#include "bvh.h"

// IGIIntegrator Method Definitions
IGIIntegrator::~IGIIntegrator() {
	
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
	
	if (dump==false) {
		return;
	}
	/*
	for (int i=0; i<virtualPaths.size(); i++) {
		for (int j=0; j<virtualPaths[i].size(); j++) {
			//printf("%s \n",virtualPaths[i][j].toString().c_str());
			dd.dump("vpl",virtualPaths[i][j].toString().c_str());
		}
	}
	for (int j=0;j<differentialRays.size();j++){
		dd.dump("rdf",differentialRays[j]->toString().c_str());
	}
	dd.dump2File(filename);
	*/
	//freeing memory arena
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

//void IGIIntegrator::setSurfaceLights(vector<VirtualSphericalLight> &vsl){
//	
//	vsls=vsl;
//	
//	int npd=vsls.size();
//	vector<Reference<Primitive> > photonDiscs;
//	photonDiscs.reserve(npd);
//	for (int i=0; i<npd; i++) {
////		Shape* shape=new PhotonDisc(new Point(),0.7);
//		//Reference<Primitive> prim=new GeometricPrimitive(shape, NULL, NULL);
//		//photonDiscs.push_back(prim) ;
//	}
//	
//
//	//printf("\n in igi vsl size is %d \n",vsls.size());
//}

void IGIIntegrator::Preprocess(const Scene *scene, const Camera *camera,
                               const Renderer *renderer) {
	
	
	
	//vsls=ps->vsls;
	
    /*if (scene->lights.size() == 0) return;
    //MC changed arena to localArena igi member variable MemoryArena arena; I have also added isFree() method to local Arena 
	if (!igiLocalArena.isFree()) {
		igiLocalArena.FreeAll();
		
	}
	//have to delete all the lights created by now
	virtualLights.clear();
	printf("\n===== vl size is %d =====\n",virtualLights.size());
    RNG rng;
	rng.Seed(1000);
    // Compute samples for emitted rays from lights
    vector<float> lightNum(nLightPaths * nLightSets);
    vector<float> lightSampPos(2 * nLightPaths * nLightSets, 0.f);
    vector<float> lightSampComp(nLightPaths * nLightSets, 0.f);
    vector<float> lightSampDir(2 * nLightPaths * nLightSets, 0.f);
    LDShuffleScrambled1D(nLightPaths, nLightSets, &lightNum[0], rng);
    LDShuffleScrambled2D(nLightPaths, nLightSets, &lightSampPos[0], rng);
    LDShuffleScrambled1D(nLightPaths, nLightSets, &lightSampComp[0], rng);
    LDShuffleScrambled2D(nLightPaths, nLightSets, &lightSampDir[0], rng);

    // Precompute information for light sampling densities
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
    for (uint32_t s = 0; s < nLightSets; ++s) {
        for (uint32_t i = 0; i < nLightPaths; ++i) {
            // Follow path _i_ from light to create virtual lights
            int sampOffset = s*nLightPaths + i;
			
            // Choose light source to trace virtual light path from
            float lightPdf;
            int ln = lightDistribution->SampleDiscrete(lightNum[sampOffset],
                                                       &lightPdf);
            Light *light = scene->lights[ln];

            // Sample ray leaving light source for virtual light path
            RayDifferential ray;
            float pdf;
            LightSample ls(lightSampPos[2*sampOffset], lightSampPos[2*sampOffset+1],
                           lightSampComp[sampOffset]);
            Normal Nl;
            Spectrum alpha = light->Sample_L(scene, ls, lightSampDir[2*sampOffset],
                                             lightSampDir[2*sampOffset+1],
                                             camera->shutterOpen, &ray, &Nl, &pdf);
			
            if (pdf == 0.f || alpha.IsBlack()) continue;
            alpha /= pdf * lightPdf;
            Intersection isect;
			
			
			
			
            while (scene->Intersect(ray, &isect) && !alpha.IsBlack()) {
                // Create virtual light and sample new ray for path
                alpha *= renderer->Transmittance(scene, RayDifferential(ray), NULL,
                                                 rng, igiLocalArena);
                Vector wo = -ray.d;
                BSDF *bsdf = isect.GetBSDF(ray, igiLocalArena);

                // Create virtual light at ray intersection point
                Spectrum contrib = alpha * bsdf->rho(wo, rng) / M_PI;
				//MC saving the bsdf in the virtual light -- for now only global radius is used 
				globalRadius=0.1;
				VirtualLight vlTemp= VirtualLight(isect.dg.p, isect.dg.nn, contrib, isect.rayEpsilon,bsdf);
				//s*i je index virtualni cesty
				virtualPaths[s*i].push_back(vlTemp);
				differentialRays.push_back(new RayDifferential(ray));
				// end MC
				
                virtualLights[s].push_back(vlTemp);

                // Sample new ray direction and update weight for virtual light path
                Vector wi;
                float pdf;
                BSDFSample bsdfSample(rng);
                Spectrum fr = bsdf->Sample_f(wo, &wi, bsdfSample, &pdf);
                if (fr.IsBlack() || pdf == 0.f)
                    break;
                Spectrum contribScale = fr * AbsDot(wi, bsdf->dgShading.nn) / pdf;

                // Possibly terminate virtual light path with Russian roulette
                float rrProb = min(1.f, contribScale.y());
                if (rng.RandomFloat() > rrProb)
                    break;
                alpha *= contribScale / rrProb;
                ray = RayDifferential(isect.dg.p, wi, ray, isect.rayEpsilon);
            }
			//MC local arena is not freed until the igi object destruction 
            //localArena.FreeAll();
        }
    }
    delete lightDistribution;
	 */
}

Spectrum IGIIntegrator::Lms(const Scene *scene, const ProgressiveRenderer *renderer,
                            const RayDifferential &ray, const Intersection &isect,
                            const Sample *sample, RNG &rng, MemoryArena &localArena) const {
    //should be called in every raycast in lss for every intersection point
    
    Spectrum ms=this->sampleVRLBruteForce(scene,renderer,ray,isect,sample,rng,localArena);
   // Spectrum ms=this->Li(scene, renderer, ray, isect, sample, rng, localArena);
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
//    rng=RNG();
//    rng.Seed(12487);
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
    // Compute indirect illumination with virtual lights
   // uint32_t lSet = min(uint32_t(sample->oneD[vlSetOffset][0] * nLightSets),
                    //    nLightSets-1);
	//MC deleted virtualLights sets
    for (uint32_t i = 0; i < this->vsls.size(); ++i) {
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
        RayDifferential connectRay(p, wi, ray, isect.rayEpsilon,
                                   sqrtf(d2) * (1.f - vl->rayEpsilon));
        Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, localArena);

        // Possibly skip virtual light shadow ray with Russian roulette
//        if (Llight.y() < rrThreshold) {
//            float continueProbability = .1f;
//            if (rng.RandomFloat() > continueProbability)
//                continue;
//            Llight /= continueProbability;
//        }

        // Add contribution from _VirtualLight_ _vl_
        if (!scene->IntersectP(connectRay))
            L += Llight;
    }
	/*
    if (ray.depth < maxSpecularDepth) {
        // Do bias compensation for bounding geometry term
        int nSamples = (ray.depth == 0) ? nGatherSamples : 1;
        for (int i = 0; i < nSamples; ++i) {
            Vector wi;
            float pdf;
            BSDFSample bsdfSample = (ray.depth == 0) ?
                BSDFSample(sample, gatherSampleOffset, i) : BSDFSample(rng);
            Spectrum f = bsdf->Sample_f(wo, &wi, bsdfSample,
                                        &pdf, BxDFType(BSDF_ALL & ~BSDF_SPECULAR));
            if (!f.IsBlack() && pdf > 0.f) {
                // Trace ray for bias compensation gather sample
                float maxDist = sqrtf(AbsDot(wi, n) / gLimit);
                RayDifferential gatherRay(p, wi, ray, isect.rayEpsilon, maxDist);
                Intersection gatherIsect;
                Spectrum Li = renderer->Li(scene, gatherRay, sample, rng, localArena,
                                           &gatherIsect);
                if (Li.IsBlack()) continue;

                // Add bias compensation ray contribution to radiance sum
                float Ggather = AbsDot(wi, n) * AbsDot(-wi, gatherIsect.dg.nn) /
                    DistanceSquared(p, gatherIsect.dg.p);
                if (Ggather - gLimit > 0.f && !isinf(Ggather)) {
                    float gs = (Ggather - gLimit) / Ggather;
                    L += f * Li * (AbsDot(wi, n) * gs / (nSamples * pdf));
                }
            }
        }
    }
	*/
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
        curVpth=vpths[i];
        //sample every vrl 5 times
        for (int s=0; s<30; ++s) {
            float d=rng.RandomFloat();
            float rayPoint=curVpth->ray.mint+(curVpth->ray.maxt-curVpth->ray.mint)*d;
            vpthPoint=curVpth->ray.o+curVpth->ray.d*rayPoint;
            Vector wi=Normalize(vpthPoint-p);
            float d2 =DistanceSquared(vpthPoint, p);
            //RayDifferential connectRay(p,wi,NULL, isect.rayEpsilon,NULL);
            RayDifferential connectRay;
            connectRay.o=p;
            connectRay.d=wi;
            connectRay.mint=isect.rayEpsilon;
           // RayDifferential connectRay(p, wi, ray, isect.rayEpsilon,NULL);
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
            //L+=curVpth->contrib*G*renderer->Transmittance(scene, RayDifferential(Ray(p, wi, 0)), NULL, rng, localArena);
           // L+=0.1*G*renderer->Transmittance(scene, RayDifferential(Ray(p, wi, 0)), NULL, rng, localArena)*renderer->Transmittance(scene, connectRay, NULL, rng, localArena);;
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
	//printf("file name is: %s",filename.c_str());
	//printf("\n params are %s \n",params.ToString().c_str());
    return new IGIIntegrator(nLightPaths, nLightSets, rrThresh,
                             maxDepth, glimit, gatherSamples,dump,filename);
	
	
}


