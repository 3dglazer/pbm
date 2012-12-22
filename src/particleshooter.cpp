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

// this method shoots light carrying paths from lights and stores these paths and when hit or scattering happens vsl is storred;
void ParticleShooter::shootParticles(const Scene * scene, Camera * camera, const Renderer *renderer, int nPaths,float radius){
	if (scene->lights.size() == 0) return;
	
	RNG rng;
	rng.Seed(seed);
	vector<float> lightNum(nPaths);
    vector<float> lightSampPos(2 * nPaths, 0.f);
    vector<float> lightSampComp(nPaths, 0.f);
    vector<float> lightSampDir(2 * nPaths, 0.f);
	
	// Precompute information for light sampling densities
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
	
	for (int currPath=0; currPath<nPaths; currPath++) {
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
		Intersection isect;	
		
		while (scene->Intersect(ray, &isect) && !alpha.IsBlack()) {
			// Create virtual light and sample new ray for path attenuate energy because of the volume transmitance
			alpha *= renderer->Transmittance(scene, RayDifferential(ray), NULL, rng, psArena);
			Vector wo = -ray.d;
			BSDF *bsdf = isect.GetBSDF(ray, psArena);
			
			// Create virtual light at ray intersection point
			Spectrum contrib = alpha * bsdf->rho(wo, rng) / M_PI;
			contrib=contrib/nPaths; //added MC
			
			//MC saving the bsdf in the virtual light -- for now only global radius is used 
			VirtualSphericalLight vslTemp= VirtualSphericalLight(isect.dg.p, isect.dg.nn, contrib, isect.rayEpsilon,radius,bsdf);
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
			ray = RayDifferential(isect.dg.p, wi, ray, isect.rayEpsilon);
		}
		//MC local arena is not freed until the igi object destruction 
		//localArena.FreeAll();
	}
	printf("\n\n======= Preprocessing %d vsls created ===========\n\n",(uint32_t)vsls.size());
	delete lightDistribution;
}



/*
// finds scatter event in homogenous media 
bool PhotonIntegrator::Transmit(const Scene * scene, RNG &rng, float t0, float t1, RayDifferential &ray, Point &interactp, Spectrum &alpha) {
    float stepsize = (rng.RandomFloat() + 0.5f)*marchSize;
    t0 += stepsize;
    Spectrum tr(1.f), step_tr, cur_tr;
	
    bool interaction = false;
	
    while (true) {
        // Stepped outside volume
        if (t0 > t1) {
            tr *= Exp(-(stepsize - (t0 - t1)) * scene->volumeRegion->sigma_t(ray(t1), Vector(), 0));
            break;
        }
		
        step_tr = Exp(-stepsize * scene->volumeRegion->sigma_t(ray(t0), Vector(), 0));
        cur_tr = tr*step_tr;
		
        // Interaction occurs
        if (rng.RandomFloat() < (1.f - cur_tr.y())) {
            // Randomly choose a point of interaction along the current interval
            float dist = (rng.RandomFloat() * stepsize);
            t0 -= dist;
            interactp = ray(t0);
            tr *= Exp(-(stepsize - dist) * scene->volumeRegion->sigma_t(ray(t0), Vector(), 0));
            interaction = true;
            break;
        }
		
        tr = cur_tr;
        
        // March with jitter
        stepsize = (rng.RandomFloat() + 0.5f)*marchSize;
        t0 += stepsize;
    }
	
    alpha *= tr;
    return interaction;
}

 
 //photon shooting in homogenous media 
 void PhotonShootingTask::Run() {
 // Declare local variables for _PhotonShootingTask_
    MemoryArena arena;
    RNG rng(31 * taskNum);
    vector<Photon> localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons;
    vector<RadiancePhoton> localRadiancePhotons;
    uint32_t totalPaths = 0;
    bool causticDone = (integrator->nCausticPhotonsWanted == 0);
    bool indirectDone = (integrator->nIndirectPhotonsWanted == 0);
    bool volumeDone = (integrator->nVolumePhotonsWanted == 0);
    PermutedHalton halton(6, rng);
    vector<Spectrum> localRpReflectances, localRpTransmittances;
    while (true) {
        // Follow photon paths for a block of samples
        const uint32_t blockSize = 4096;
        for (uint32_t i = 0; i < blockSize; ++i) {
            float u[6];
            halton.Sample(++totalPaths, u);
            // Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];
			
            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);
            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
            if (!alpha.IsBlack()) {
                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);
                bool specularPath = true;
                Intersection photonIsect;
                int nIntersections = 0;
                bool absorbed = false;
                int nScatters = 0;
                
                while (scene->Intersect(photonRay, &photonIsect)) {
                    float t0, t1;
                    bool surfacehit = true;
					
                    // Handle photon/volume intersection
                    while (scene->volumeRegion->IntersectP(photonRay, &t0, &t1)) {
                        Point interactp;
                        
                        if (integrator->Transmit(scene, rng, t0, t1, photonRay, interactp, alpha)) {
                            Spectrum albedo = scene->volumeRegion->sigma_s(interactp, Vector(), 0.f)/scene->volumeRegion->sigma_t(interactp, Vector(), 0.f);
							
                            // Store photon only if indirect light
                            if (!volumeDone && (nScatters > 0 || nIntersections > 0)) {
                                Photon photon(interactp, alpha, photonRay.d);
                                localVolumePhotons.push_back(photon);
                                //photonsfile << "AttributeBegin\n";
                                //photonsfile << "Material \"matte\" \"color Kd\" [0 1 0]\n";
                                //photonsfile << "Translate " << interactp.x << " " << interactp.y << " " << interactp.z << "\n";
                                //photonsfile << "Shape \"sphere\" \"float radius\" [.005]\n";
                                //photonsfile << "AttributeEnd\n\n";
                            }
							
                            // Scatter
                            if (rng.RandomFloat() < albedo.y()) {
                                ++nScatters;
                                // Sample new direction
                                Vector newdir = SampleHG(photonRay.d, scene->volumeRegion->g, rng.RandomFloat(), rng.RandomFloat());
                                alpha *= scene->volumeRegion->p(interactp, photonRay.d, newdir, 0.f);
                                alpha *= albedo;
								
                                photonRay = RayDifferential(interactp, newdir, 0.f);
                                scene->Intersect(photonRay, &photonIsect);
                            } else { // Absorb
                                surfacehit = false;
                                break;
                            }
                        } else {
                            photonRay = RayDifferential(photonRay(t1), photonRay.d, 0.f);
                            break;
                        }
                    }
					
                    if (!surfacehit) break;
					
                    ++nIntersections;
					
                    // Handle photon/surface intersection
                    alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
					
                    BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
                    BxDFType specularType = BxDFType(BSDF_REFLECTION |
													 BSDF_TRANSMISSION | BSDF_SPECULAR);
                    bool hasNonSpecular = (photonBSDF->NumComponents() >
										   photonBSDF->NumComponents(specularType));
                    Vector wo = -photonRay.d;
                    if (hasNonSpecular) {
                        
                        // Deposit photon at surface
                        Photon photon(photonIsect.dg.p, alpha, wo);
                        bool depositedPhoton = false;
                        if (specularPath && nIntersections > 1) {
                            if (!causticDone) {
                                PBRT_PHOTON_MAP_DEPOSITED_CAUSTIC_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true;
                                localCausticPhotons.push_back(photon);
                            }
                        }
                        else {
                            // Deposit either direct or indirect photon
                            // stop depositing direct photons once indirectDone is true; don't
                            // want to waste memory storing too many if we're going a long time
                            // trying to get enough caustic photons desposited.
                            if (nIntersections == 1 && !indirectDone && integrator->finalGather) {
                                PBRT_PHOTON_MAP_DEPOSITED_DIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true;
                                localDirectPhotons.push_back(photon);
                            }
                            else if (nIntersections > 1 && !indirectDone) {
                                PBRT_PHOTON_MAP_DEPOSITED_INDIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                                depositedPhoton = true;
                                localIndirectPhotons.push_back(photon);
                            }
                        }
						
                        // Possibly create radiance photon at photon intersection point
                        if (depositedPhoton && integrator->finalGather &&
							rng.RandomFloat() < .125f) {
                            Normal n = photonIsect.dg.nn;
                            n = Faceforward(n, -photonRay.d);
                            localRadiancePhotons.push_back(RadiancePhoton(photonIsect.dg.p, n));
                            Spectrum rho_r = photonBSDF->rho(rng, BSDF_ALL_REFLECTION);
                            localRpReflectances.push_back(rho_r);
                            Spectrum rho_t = photonBSDF->rho(rng, BSDF_ALL_TRANSMISSION);
                            localRpTransmittances.push_back(rho_t);
                        }
                    }
					
                    if (nIntersections >= integrator->maxPhotonDepth) break;
					
                    // Sample new photon ray direction
                    Vector wi;
                    float pdf;
                    BxDFType flags;
                    Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
													   &pdf, BSDF_ALL, &flags);
                    if (fr.IsBlack() || pdf == 0.f) break;
                    Spectrum anew = alpha * fr *
					AbsDot(wi, photonBSDF->dgShading.nn) / pdf;
					
                    // Possibly terminate photon path with Russian roulette
                    float continueProb = min(1.f, anew.y() / alpha.y());
                    if (rng.RandomFloat() > continueProb)
                        break;
                    alpha = anew / continueProb;
                    specularPath &= ((flags & BSDF_SPECULAR) != 0);
                    
                    if (indirectDone && !specularPath) break;
                    photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
                                                photonIsect.rayEpsilon);
                    
                }
                PBRT_PHOTON_MAP_FINISHED_RAY_PATH(&photonRay, &alpha);
            }
            arena.FreeAll();
        }
		
        // Merge local photon data with data in _PhotonIntegrator_
        { MutexLock lock(mutex);
			
			// Give up if we're not storing enough photons
			if (abortTasks)
				return;
			if (nshot > 500000 &&
				(unsuccessful(integrator->nCausticPhotonsWanted,
							  causticPhotons.size(), blockSize) ||
				 unsuccessful(integrator->nIndirectPhotonsWanted,
							  indirectPhotons.size(), blockSize) || 
				 unsuccessful(integrator->nVolumePhotonsWanted,
							  volumePhotons.size(), blockSize))) {
					 Error("Unable to store enough photons.  Giving up.\n");
					 causticPhotons.erase(causticPhotons.begin(), causticPhotons.end());
					 indirectPhotons.erase(indirectPhotons.begin(), indirectPhotons.end());
					 volumePhotons.erase(volumePhotons.begin(), volumePhotons.end());
					 radiancePhotons.erase(radiancePhotons.begin(), radiancePhotons.end());
					 abortTasks = true;
					 return;
				 }
			progress.Update(localIndirectPhotons.size() + localCausticPhotons.size() + localVolumePhotons.size());
			nshot += blockSize;
			
			// Merge indirect photons into shared array
			if (!indirectDone) {
				integrator->nIndirectPaths += blockSize;
				for (uint32_t i = 0; i < localIndirectPhotons.size(); ++i)
					indirectPhotons.push_back(localIndirectPhotons[i]);
				localIndirectPhotons.erase(localIndirectPhotons.begin(),
										   localIndirectPhotons.end());
				if (indirectPhotons.size() >= integrator->nIndirectPhotonsWanted)
					indirectDone = true;
				nDirectPaths += blockSize;
				for (uint32_t i = 0; i < localDirectPhotons.size(); ++i)
					directPhotons.push_back(localDirectPhotons[i]);
				localDirectPhotons.erase(localDirectPhotons.begin(),
										 localDirectPhotons.end());
			}
			
			// Merge direct, caustic, and radiance photons into shared array
			if (!causticDone) {
				integrator->nCausticPaths += blockSize;
				for (uint32_t i = 0; i < localCausticPhotons.size(); ++i)
					causticPhotons.push_back(localCausticPhotons[i]);
				localCausticPhotons.erase(localCausticPhotons.begin(), localCausticPhotons.end());
				if (causticPhotons.size() >= integrator->nCausticPhotonsWanted)
					causticDone = true;
			}
			
			// Merge volume photons
			if (!volumeDone) {
				integrator->nVolumePaths += blockSize;
				for (uint32_t i = 0; i < localVolumePhotons.size(); ++i)
					volumePhotons.push_back(localVolumePhotons[i]);
				localVolumePhotons.erase(localVolumePhotons.begin(), localVolumePhotons.end());
				if (volumePhotons.size() >= integrator->nVolumePhotonsWanted)
					volumeDone = true;
			}
			
			//printf("%d\n", indirectPhotons.size());
			
			for (uint32_t i = 0; i < localRadiancePhotons.size(); ++i)
				radiancePhotons.push_back(localRadiancePhotons[i]);
			localRadiancePhotons.erase(localRadiancePhotons.begin(), localRadiancePhotons.end());
			for (uint32_t i = 0; i < localRpReflectances.size(); ++i)
				rpReflectances.push_back(localRpReflectances[i]);
			localRpReflectances.erase(localRpReflectances.begin(), localRpReflectances.end());
			for (uint32_t i = 0; i < localRpTransmittances.size(); ++i)
				rpTransmittances.push_back(localRpTransmittances[i]);
			localRpTransmittances.erase(localRpTransmittances.begin(), localRpTransmittances.end());
        }
		
        // Exit task if enough photons have been found
        if (indirectDone && causticDone && volumeDone)
            break;
    }
}


*/