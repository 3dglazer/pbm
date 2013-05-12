
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


// main/pbrt.cpp*
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"
#include "samplingmethods.h"
// main program
int main(int argc, char *argv[]) {
    
    Point vrl0=Point(0,0,0);
    Point vrlend=Point(1,0,0);
    
    Point p=Point(0.5,1.,0.);
    float sParam;
    float h=p2Ray(p, vrl0, vrlend, sParam);
    float s=0;
    printf("\nh=%f; sParam=%f;\n",h,sParam);
    float rayPoint=invEqAngCDF(s, h, 0.-sParam, 1.-sParam);
    float u0=0.-sParam;
    float u1=1.-sParam;
    float val0=atanf(u0/h);
    float val1=atanf(u1/h);
    float lrp1=Lerp(s, val0, val1);
    float lrp0=Lerp(val0, val1, s);
    
    float newRayPoint=h*tan(s);
    
    printf("\n newRayPoint=%f",newRayPoint);
    
    
    Options options;
    vector<string> filenames;
    // Process command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--ncores")) options.nCores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--outfile")) options.imageFile = argv[++i];
        else if (!strcmp(argv[i], "--quick")) options.quickRender = true;
        else if (!strcmp(argv[i], "--quiet")) options.quiet = true;
        else if (!strcmp(argv[i], "--verbose")) options.verbose = true;
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: pbrt [--ncores n] [--outfile filename] [--quick] [--quiet] "
                   "[--verbose] [--help] <filename.pbrt> ...\n");
            return 0;
        }
        else filenames.push_back(argv[i]);
    }

    // Print welcome banner
    if (!options.quiet) {
        printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
               PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
        printf("Copyright (c)1998-2010 Matt Pharr and Greg Humphreys.\n");
        printf("The source code to pbrt (but *not* the book contents) is covered by the GNU GPL.\n");
        printf("See the file COPYING.txt for the conditions of the license.\n");
        fflush(stdout);
    }
    pbrtInit(options);
    // Process scene description
    PBRT_STARTED_PARSING();
    if (filenames.size() == 0) {
        // Parse scene from standard input
        ParseFile("-");
    } else {
        // Parse scene from input files
        for (u_int i = 0; i < filenames.size(); i++)
            if (!ParseFile(filenames[i]))
                Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
    }
    pbrtCleanup();
    return 0;
}


