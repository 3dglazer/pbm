//
//  main.cpp
//  PHOTONVIZ
//
//  Created by Jaromir Herskovic on 06.05.13.
//  Copyright (c) 2013 ZDENEK GLAZER. All rights reserved.
//

#include <iostream>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "PCloud.h"
#include <stdlib.h>
#include "Parser.h"
#include "OGLRenderer.h"
#include "mystructs.h"
int main(int argc, char * argv[])
{
    glutInit( &argc, argv );
    PCloud pc;
    std::vector<float*> points;
    vector<volumePath *> paths;

    //pc.generateGaussianDistribution(points, 100, 3);
    Parser p;
    
    p.getPoints(points, "/Volumes/DISK2/developer/pbm/photonVIZ/pointCloud.txt");
    p.getLines(paths,"/Volumes/DISK2/developer/pbm/photonVIZ/volumePaths.txt");
//    for (int i=0; i<paths.size(); i++) {
//        printf("\nVpath:");
//        volumePath *vp=paths[i];
//        for (int j=0; j<vp->points.size(); j++) {
//            float * psr=vp->points[j];
//            printf("\n%f,%f,%f,%f\n",psr[0],psr[1],psr[2],psr[3]);
//        }
//    }
    
    OGLRenderer renderer(10.);
    if (!paths.empty()) {
        printf("PATHS ARE EMPTY!\n");
        renderer.setVolumeLines(paths);
    }
    if (!points.empty()) {
        printf("POINTS ARE EMPTY!\n");
        renderer.setPoints(points);
    }
    renderer.resetSceneCamera();
    //renderer.setupGLUTCallbacks();
    renderer.render();
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

