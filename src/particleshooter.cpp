/*
 *  particleshooter.cpp
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/18/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */

#include "particleshooter.h"
#include "montecarlo.h"
#include "sampler.h"
#include "intersection.h"
#include <algorithm>

// this method shoots light carrying paths from lights and stores these paths and when hit or scattering happens vsl is storred;
void ParticleShooter::shootParticles(const Scene * scene, Camera * camera, const ProgressiveRenderer *renderer, int nPaths,float radius){
    bool hasVolumes=true;
    if (scene->volumeRegion) {
        printf("\n Photon caches are being destributed over volume and surfices.\n");
    }else{
        hasVolumes=false;
        printf("\n Photon caches are being destributed only over surfices no volumes found!!\n");   
    }
	if (scene->lights.size() == 0) return;
	
	RNG rng;
	rng.Seed(seed);
    int nTimesTries=3; //how many times more samples will be tried to be shoot to get nPaths
	vector<float> lightNum(nPaths*nTimesTries);
    vector<float> lightSampPos(2 * nTimesTries * nPaths, 0.f);
    vector<float> lightSampComp(nPaths * nTimesTries, 0.f);
    vector<float> lightSampDir(2 * nTimesTries *nPaths, 0.f);
    LDShuffleScrambled1D(nPaths*nTimesTries, 1, &lightNum[0], rng);
    LDShuffleScrambled2D(nPaths*nTimesTries, 1, &lightSampPos[0], rng);
    LDShuffleScrambled1D(nPaths*nTimesTries, 1, &lightSampComp[0], rng);
    LDShuffleScrambled2D(nPaths*nTimesTries, 1, &lightSampDir[0], rng);
	// Precompute information for light sampling densities
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
	
	for (int currPath=0; vsls.size()<=nPaths || (currPath<nPaths*nTimesTries); currPath++) {
		// Choose light source to trace virtual light path from
		float lightPdf;
		int ln = lightDistribution->SampleDiscrete(lightNum[currPath],
												   &lightPdf);
		Light *light = scene->lights[ln];
		
		// Sample ray leaving light source for virtual light path
		RayDifferential ray;
		float pdf;
		LightSample ls(lightSampPos[2*currPath], lightSampPos[2*currPath+1],
					   lightSampComp[currPath]);
		Normal Nl;
		Spectrum alpha = light->Sample_L(scene, ls, lightSampDir[2*currPath],
										 lightSampDir[2*currPath+1],
										 camera->shutterOpen, &ray, &Nl, &pdf);
		
		if (pdf == 0.f || alpha.IsBlack()) continue;
		alpha /= pdf * lightPdf;
        alpha /=nPaths;
		Intersection isect;	
		
        //check whether intersect any scene primitive, including volumeRegions
		while (scene->Intersect(ray, &isect)&&!alpha.IsBlack()) {
            float isectDist=DistanceSquared(ray.o, isect.dg.p);
            //perform Woodcock tracking
            float t0,t1;
            Spectrum tau(0.);
            float d;
            //perform intersection of volume agregate, gives us the t0 from, t1 to distance parametres 
            if (hasVolumes && scene->volumeRegion->IntersectP(ray, &t0, &t1)) {
                //stop criteria in freeFlight are ray.mint and ray. maxt
                if (isectDist<t0*t0)continue;
                if (isectDist<t1*t1)t1=sqrtf(isectDist);
                ray.mint=t0;
                ray.maxt=t1;
                //create a VRL the ray vector has to be normalizet otherwise wont work!!!!
                ray.d=Normalize(ray.d);
                VolumePath* currPath=new VolumePath(ray);
                currPath->contrib=alpha; //sets the energy to the path
                volumePaths.push_back(currPath);
                //maybe add some constraining criteria like maximum scattering events in the original VRL paper is 16
                Spectrum maxtau;
                float dmax=-INFINITY;
                for (int evnts=0; evnts < maxScattering; ++evnts) {
                    d=renderer->freeFlight(scene, ray, maxtau, rng); //the tau could be used for multiple scattering
                    if (d==-1.){
                        continue; //scattering did not happened tau should be valid
                    }else{
                        if (d>dmax) {
                            dmax=d;
                            tau=maxtau;
                        }
                        //ray.mint=d; //save event distance for next tracking
                        //add the transmittance and distance to the VRL
                        currPath->dists.push_back(d);
                        //phase function doesn't have to be saved can be querried directly from scene->volumeRegion->p(....) returns probability float.
                    }
                }// end of SCATTERING EVENTS
                //if the tracking ended before reaching the end point of media we have to add the optical thickness of the rest of the media to it, to correctly attenuate the vsls.
                if (dmax<t1) {
                    alpha*=renderer->Transmittance(scene, ray, NULL, rng, psArena); //attnuation to the vsl
                    //tau+=(t1-d)*scene->volumeRegion->sigma_t(ray.o+ray.d*t1, ray.d, ray.time);
                }
                std::sort((currPath->dists.begin()), (currPath->dists.end()));
                printf("\n");
                for (int i=0; i<currPath->dists.size(); ++i) {
                    printf("%f, ",currPath->dists[i]);
                }
            }
            //VRL has been created only if volume interaction occured not it is time to create VPL
            //MC commented
			//alpha *= renderer->Transmittance(scene, RayDifferential(ray), NULL, rng, psArena);
            alpha*=Exp(-tau); //transmittance from optical thickness
			Vector wo = -ray.d;
			BSDF *bsdf = isect.GetBSDF(ray, psArena);
			
			// Create virtual light at ray intersection point
			//Spectrum contrib = alpha * bsdf->rho(wo, rng) / M_PI; //Lambertian contribution
			Spectrum contrib=alpha; //added MC not very clever all paths might not be created in all, no nead to devide it by nPahts, because it has been scaled in the beginning by pdf
			
			//MC saving the bsdf in the virtual light -- for now only global radius is used 
			VirtualSphericalLight *vslTemp= new VirtualSphericalLight(isect.dg.p,wo,isect.dg.nn, contrib, isect.rayEpsilon,radius,bsdf);
			//s*i je index virtualni cesty
			//virtualPaths[s*i].push_back(vlTemp);
			//differentialRays.push_back(new RayDifferential(ray));
			// end MC
			
			vsls.push_back(vslTemp);
			
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
            //alpha *=contribScale;
			ray = RayDifferential(isect.dg.p, wi, ray, isect.rayEpsilon);
		}
		//MC local arena is not freed until the igi object destruction 
		//localArena.FreeAll();
	}
	printf("\n\n======= Preprocessing %d vsls created ===========\n\n",(uint32_t)vsls.size());
	delete lightDistribution;
}
/*
 //old woodcock tracking
//stop criteria in freeFlight are ray.mint and ray. maxt
if (isectDist<t0*t0)continue;
if (isectDist<t1*t1)t1=sqrtf(isectDist);
ray.mint=t0;
ray.maxt=t1;
//create a VRL the ray vector has to be normalizet otherwise wont work!!!!
ray.d=Normalize(ray.d);
VolumePath* currPath=new VolumePath(ray);
currPath->contrib=alpha; //sets the energy to the path
volumePaths.push_back(currPath);
//maybe add some constraining criteria like maximum scattering events in the original VRL paper is 16
for (int evnts=0; evnts < maxScattering; ++evnts) {
    d=renderer->freeFlight(scene, ray, tau, rng); //the tau could be used for multiple scattering
    if (d==-1.){
        break; //scattering did not happened tau should be valid
    }else{
        ray.mint=d; //save event distance for next tracking
        //add the transmittance and distance to the VRL
        currPath->dists.push_back(d);
        //phase function doesn't have to be saved can be querried directly from scene->volumeRegion->p(....) returns probability float.
    }
    //possibly terminate tracking if the transmittance is really small
    if (tau.returnOne() < 1e-3) {
        const float continueProb = .5f;
        if (rng.RandomFloat() > continueProb) break;
    }
}// end of SCATTERING EVENTS
//if the tracking ended before reaching the end point of media we have to add the optical thickness of the rest of the media to it, to correctly attenuate the vsls.
if (d<t1) {
    tau+=(t1-d)*scene->volumeRegion->sigma_t(ray.o+ray.d*t1, ray.d, ray.time);
}
*/