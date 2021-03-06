
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_RENDERER_H
#define PBRT_CORE_RENDERER_H

// core/renderer.h*
#include "pbrt.h"

// Renderer Declarations
class Renderer {
public:
    // Renderer Interface
    virtual ~Renderer();
    virtual void Render(const Scene *scene) = 0;
    
    virtual Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const = 0;
    
    virtual Spectrum Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample,
        RNG &rng, MemoryArena &arena) const = 0;
    


};

class ProgressiveRenderer: public Renderer{
    public:
        virtual Spectrum Lms(const Scene *scene, const RayDifferential &ray, const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect, Spectrum *T) const =0;
        virtual Spectrum Lss(const Scene *scene,
                     const RayDifferential &ray, const Sample *sample, RNG &rng,
                     MemoryArena &arena, Intersection *isect, Spectrum *T) const = 0;
        virtual Spectrum Lsm(const Scene *scene,
                     const RayDifferential &ray, const Sample *sample, RNG &rng,
                     MemoryArena &arena, Intersection *isect, Spectrum *T) const = 0;
        virtual Spectrum Lmm(const Scene *scene,
                     const RayDifferential &ray, const Sample *sample, RNG &rng,
                     MemoryArena &arena, Intersection *isect, Spectrum *T) const = 0;
        virtual float freeFlight(const Scene *scene, const Ray &r,Spectrum& tau,const RNG &rng) const = 0;
};

#endif // PBRT_CORE_RENDERER_H
