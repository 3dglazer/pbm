
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


// renderers/samplerrenderer.cpp*
#include "stdafx.h"
#include "renderers/progressiveRenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"

static uint32_t hash(char *key, uint32_t len)
{
    uint32_t   hash, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

void PGRendererTask::repairSample(Spectrum &s){
    // Issue warning if unexpected radiance value returned
    if (s.HasNaNs()) {
//        Error("Not-a-number radiance value returned "
//              "for image sample.  Setting to black.");
        s = Spectrum(0.f);
    }
    else if (s.y() < -1e-5) {
//        Error("Negative luminance value, %f, returned"
//              "for image sample.  Setting to black.", s.y());
        s = Spectrum(0.f);
    }
    else if (isinf(s.y())) {
//        Error("Infinite luminance value returned"
//              "for image sample.  Setting to black.");
        s = Spectrum(0.f);
    }
}
// PGRendererTask Definitions
void PGRendererTask::Run() {

    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }
	
    // Declare local variables used for rendering loop
    MemoryArena arena;
	//MC added seed value;
    RNG rng(taskNum*seed);
	
    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples]; //was contribution from Li so final image now too.
    Spectrum *Ts = new Spectrum[maxSamples];
    //MC added samples for different film plates
    Spectrum *Lss = new Spectrum[maxSamples];
    Spectrum *Lsm = new Spectrum[maxSamples];
    Spectrum *Lmm = new Spectrum[maxSamples];
    Spectrum *Lms = new Spectrum[maxSamples];
    //end of MC
    Intersection *isects = new Intersection[maxSamples];
	float rgbdef[3]={0.,0.,0.};
    Spectrum zero=Spectrum::FromRGB(rgbdef);
    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);
			
            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            if (visualizeObjectIds) {
                if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
                    // random shading based on shape id...
                    uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                    uint32_t h = hash((char *)ids, sizeof(ids));
                    float rgb[3] = { (h & 0xff), (h >> 8) & 0xff, (h >> 16) & 0xff };
                    Ls[i] = Spectrum::FromRGB(rgb);
                    Ls[i] /= 255.f;
                }
                else
                    //if the ray doesnt intersect geometry
//                    for (uint32_t j = 0; j < scene->lights.size(); ++j){
//                        Li += scene->lights[j]->Le(rays[i]);
//                    }
//                    Lss[i]=Li;
                    Ls[i] = 0.f;
                
            }
            else {
				if (rayWeight > 0.f){
                    
                    Ls[i]=zero;
                    Lss[i]=zero;
                    Lsm[i]=zero;
                    Lms[i]=zero;
                    Lmm[i]=zero;
                    Intersection *isect=&isects[i];
                    RayDifferential currRay=rays[i];
                    currRay.maxt=10000.f;
                    Sample *currSample=&samples[i];
                    Spectrum tempT=1.;
                    Spectrum Li;
                    //printf("weightGreater than zero\n");
//                    Spectrum rd=RGBSpectrum();
//                    float* rcl=new float[3];
//                    rcl[0]=1.0;
//                    rcl[1]=0.;
//                    rcl[2]=0.;
//                    rd.FromRGB(&rcl[0]);
                    //printf("%f,%f,%f",rd);
                    //scene->Intersect(currRay, isect);
                    //Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng, arena, &isects[i], &Ts[i]);
                    if (scene->Intersect(rays[i], &isects[i])) {
                        //printf("intersection happened");
                        //Ls[i] = rayWeight * surfaceIntegrator->Li(scene, renderer, rays[i], isects[i], currSample,rng, arena);
                        Lss[i] = rayWeight * surfaceIntegrator->Li(scene, renderer, currRay, *isect, currSample,rng, arena);
                        //Lss[i]=rd;
                        Lms[i] = rayWeight * surfaceIntegrator->Lms(scene, renderer, currRay, *isect, currSample,rng, arena);
                    }else{
                        // Handle ray that doesn't intersect any geometry

                    }
                    Lmm[i]= rayWeight * volumeIntegrator->Lmm(scene, renderer, RayDifferential(currRay), currSample, rng, &Ts[i], arena);
                    Spectrum tr=Ts[i];
                    //Lmm[i]= rayWeight * volumeIntegrator->Li(scene, renderer, currRay, currSample, rng, &Ts[i], arena);
                    Lsm[i]= rayWeight * volumeIntegrator->Lsm(scene, renderer, RayDifferential(currRay), currSample, rng, &tempT, arena);
                    //Ts[i]=tempT;
                    Ls[i]= (Lss[i]+Lms[i])*(Ts[i])+ Lmm[i] + Lsm[i];
                }else {
					Ls[i] = 0.f;
                    //MC
                    Lmm[i]= 0.f;
                    Lsm[i]= 0.f;
                    Lss[i]= 0.f;
                    Lms[i]= 0.f;
                    //end fo MC
					Ts[i] = 1.f;
				}
                //should repair samples here
                repairSample(Ls[i]);
                repairSample(Lmm[i]);
                repairSample(Lms[i]);
                repairSample(Lsm[i]);
                repairSample(Lss[i]);
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }
		
        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i]);
				//MC added progCamera
				progCamera->film->AddSample(samples[i], Ls[i]);
				mm->AddSample(samples[i], Lmm[i]);
                ms->AddSample(samples[i], Lms[i]);
                sm->AddSample(samples[i], Lsm[i]);
                ss->AddSample(samples[i], Lss[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }
		
        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }
	
    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
								sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Lmm;
    delete[] Lsm;
    delete[] Lss;
    delete[] Lms;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}


// SamplerRenderer Method Definitions
PGRenderer::PGRenderer(Sampler *s, Camera *c,Camera *prc,Film* surface2surface,Film* surface2media,Film* media2surface, Film* media2media, ProgressiveSurfaceIntegrator *si,
										 ProgressiveVolumeIntegrator *vi, bool visIds,int nIterations,int nps,float rad,int rndSeed) {
    sampler = s;
    camera = c;
	progCamera=prc;
	//films
	ss=surface2surface;
	ms=media2surface;
	mm=media2media;
	sm=surface2media;
	//end of films
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
	nIter=nIterations;
	seed=rndSeed;
//MC for particle shooter 
	nParticles=nps;
	radius=rad;
}


PGRenderer::~PGRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}

//MC
float PGRenderer::freeFlight(const Scene *scene, const Ray &r,Spectrum& tau,const RNG &rng) const{
    return volumeIntegrator->freeFlight(scene, r, tau, rng);
}


//MC this method does every iteration rendering
void PGRenderer::renderIter(int currentIter,const Scene *scene, Sample *sample){
	// Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
	//MC added particle shooter, which is preprocessing stage for photon or vpls,vsl,vrl shooting
	ParticleShooter * particleShooter=new ParticleShooter(seed*((currentIter+1)*1372568),8);
	printf("Preprocessing stage shooting particles nPaths %d, %f",nParticles,radius);
	particleShooter->shootParticles(scene, camera, this,nParticles,radius);
    particleShooter->dumpVSLS("/Volumes/DISK2/developer/pbm/photonVIZ/pointCloud.txt");
     particleShooter->dumpVpths("/Volumes/DISK2/developer/pbm/photonVIZ/volumePaths.txt");
    printf("after particle shooting");
	surfaceIntegrator->setSurfaceLights(particleShooter->vsls);
    surfaceIntegrator->setVolumeLights(particleShooter->volumePaths);
    surfaceIntegrator->Preprocess(scene, camera, this);
    //-------- TODO ----------
    //set vrls made in particle shooting preprocess
    volumeIntegrator->setSurfaceLights(particleShooter->vsls);
    volumeIntegrator->setVolumeLights(particleShooter->volumePaths);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();

	
    // Create and launch _SamplerRendererTask_s for rendering image
	
    // Compute number of _SamplerRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
    nTasks = RoundUpPow2(nTasks);
	
	char num[16]; // string which will contain the number
	sprintf ( num, "%d", currentIter );
	string rendr= "Rendering_frameNumber"+string(num);
    ProgressReporter reporter(nTasks,rendr);
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i)
		//MC added progressive camera
        renderTasks.push_back(new PGRendererTask(scene, this,this->volumeIntegrator,this->surfaceIntegrator, camera,progCamera,ss,sm,ms,mm,
														  reporter, sampler, sample, 
														  visualizeObjectIds, 
														  nTasks-1-i, nTasks,seed*(currentIter+1)));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();
    
	//MC adding frame number after name 
    camera->film->WriteIterImage(currentIter);
//    mm->WriteIterImage(currentIter);
//    ms->WriteIterImage(currentIter);
//    ss->WriteIterImage(currentIter);
//    sm->WriteIterImage(currentIter);
	delete particleShooter;//have to be deleted last, because of the memory arena which holds brdf information from the virtual light shooting
	// have to 
}

void PGRenderer::Render(const Scene *scene) {

    PBRT_FINISHED_PARSING();
	// Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);

	for (int n=0; n<nIter; n++) {
		renderIter(n, scene,sample);
		//delete all the information by now
		//camera->film->resetPixels(); //might cause problem originaly was not commented and displayed individual iterations
	}

    //MC write all the different transport averaged images to separate files
	mm->WriteImage();
    ms->WriteImage();
    ss->WriteImage();
    sm->WriteImage();
    //write averaged image
	progCamera->film->WriteImage();
	// Clean up after rendering and store final image
    delete sample;
	PBRT_FINISHED_RENDERING();
	
}


// ve volume integratoru
//Spectrum ::Lmm
//Spectrum ::Lsm
// ve vrl integratoru
//Spectrum ::Lms
//Spectrum ::Lss


Spectrum PGRenderer::Li(const Scene *scene,
							 const RayDifferential &ray, const Sample *sample, RNG &rng,
							 MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect)){
    Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena);
    }else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
			Li += scene->lights[i]->Le(ray);
    }
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    return *T * Li + Lvi;
    
}



Spectrum PGRenderer::Transmittance(const Scene *scene,
										const RayDifferential &ray, const Sample *sample, RNG &rng,
										MemoryArena &arena) const {

    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}
//These methods are for external use only
Spectrum PGRenderer::Lms(const Scene *scene,
                         const RayDifferential &ray, const Sample *sample, RNG &rng,
                         MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect)){
        Li = surfaceIntegrator->Lms(scene, this, ray, *isect, sample,
                                    rng, arena);
    }
    return *T * Li ;
}

Spectrum PGRenderer::Lss(const Scene *scene,
                         const RayDifferential &ray, const Sample *sample, RNG &rng,
                         MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect)){
        Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                    rng, arena);
    }else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
			Li += scene->lights[i]->Le(ray);
    }
    //return *T * Li ;
    return Li;
    
}

Spectrum PGRenderer::Lsm(const Scene *scene,
                         const RayDifferential &ray, const Sample *sample, RNG &rng,
                         MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect)){
        Li = volumeIntegrator->Lsm(scene, this, ray, sample, rng,
                                      T, arena);
    }
    return  Li ;
    
}

Spectrum PGRenderer::Lmm(const Scene *scene,
                         const RayDifferential &ray, const Sample *sample, RNG &rng,
                         MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect)){
        Li = volumeIntegrator->Lmm(scene, this, ray, sample, rng,
                                      T, arena);
    }
    return  Li ;
    
}

