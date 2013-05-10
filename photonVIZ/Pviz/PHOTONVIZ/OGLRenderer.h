/*
 *  OpenGLRenderer.h
 *  shadowMapping
 *
 *  Created by System Administrator on 2/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
//Header lock, prevents multiple include build errors
#ifndef OGL_RENDERER_H
#define OGL_RENDERER_H

#ifdef _WIN32
#include "windows.h"
#else
//#include "GLUT/glut.h"
#endif
#include <cstdlib>
#include <stdio.h>
#include "BBox.h"
//#include "KDTree.h"
#include <vector>
#include <stack>
#include "mystructs.h"
using namespace std;

/// this class is used for kdtree visualization
class OGLRenderer {
public:
    //constructor
	OGLRenderer(double sceneDiametr);
	
	//destructor
	~OGLRenderer(){};
    void resetSceneCamera();
    void setPoints(vector<float*> points);
    void CameraRefresh();
    void setVolumeLines(vector<volumePath*> &paths){
        this->_volPaths=paths;
//        for (int i=0; i<paths.size(); i++) {
//            printf("\nVpath:");
//            volumePath *vp=this->_volPaths[i];
//            for (int j=0; j<vp->points.size(); j++) {
//                float * psr=vp->points[j];
//                printf("\n%f,%f,%f,%f\n",psr[0],psr[1],psr[2],psr[3]);
//            }
//        }
    };
private:
	int windowWidth;
	int windowHeight;
	
	// mouse controls
	int mouse_old_x, mouse_old_y;
	int mouse_buttons;
	float translate_z;
    float defaultScale;
	float rotate_x , rotate_y;
	vector<BBox*> _bBoxes;
	vector<float*> _origPoints;
    vector<volumePath*> _volPaths;
	//Node* _root;
	void drawBox( float min[3], float max[3]);
	void drawLight();
	void drawBBoxes();
	void drawPoints();
	//void debugTree(Node& currNode, BBox& currBBox);
	int _depth;
	int _maxDepthToRender;
	public:
	//GLUT CAN'T BE CALLED DIRECTLY FROM C++ CLASS FUNCTIONS, HACK USED. MORE INFRO ON http://stackoverflow.com/questions/3589422/using-opengl-glutdisplayfunc-within-class
	void setupGLUTCallbacks();
	void drawScene();
	void reloadPoints();
	void setBBoxes(vector<BBox*> bBoxes);
	//void setRoot(Node& root);
	void handleMouse(int button, int state, int x, int y);
	void handleKeys(unsigned char key, int x, int y);
	void cameraMotion(int x, int y);
	void idleFunc();
	void render();
    void drawLines();

	

};
#endif