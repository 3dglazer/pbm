/*
 *  OpenGLRenderer.cpp
 *  shadowMapping
 *
 *  Created by System Administrator on 2/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#define TREECOLOR 0.4f,0.4f,0.4f
#define CURRTREECOLOR 0.6f,1.0f,0.4f
#define	VISITEDTREECOLOR  0.5f,0.5f,1.0f
#define POINTCOLOR 0.2f,0.5f,0.2f
#define CURRPOINTCOLOR 1.0f,0.5f,0.5f
#define VISITEDPOINTCOLOR  0.5f,1.0f,1.0f
#define IMAGE_WIDTH  1024
#define IMAGE_HEIGHT 1024

#include "OGLRenderer.h"
#include <GLUT/glut.h>
#include <cmath>
#include <math.h>

float randomF(){
	float scale=RAND_MAX+1.;
	float base=rand()/scale;
	float fine=rand()/scale;
	return base+fine/scale;
}

//============== GLUT HACK ==================
//GLUT CAN'T BE CALLED DIRECTLY FROM C++ CLASS FUNCTIONS, HACK USED.
//MORE INFRO ON http://stackoverflow.com/questions/3589422/using-opengl-glutdisplayfunc-within-class
OGLRenderer* currentInstance;
extern "C"
void drawCallback()
{
	currentInstance->drawScene();
}
// static function for draw callback for GLUT and OpenGL
extern "C"
void idleCallback(){
	currentInstance->idleFunc();
}
extern "C"
void mouseCallback(int button, int state, int x, int y){
	currentInstance->handleMouse(button,state,x,y);
}
extern "C"
void motionCallback(int x , int y){
	currentInstance->cameraMotion(x,y);
}
extern "C"
void keyCallback(unsigned char key, int x, int y){
	currentInstance->handleKeys(key,x,y);	
}
//=========== END OF GLUT HACK ===============
void OGLRenderer::resetSceneCamera(){
	double diam=20.;
	double zNear = 0.00001;
	double zFar = zNear + diam*5;
	GLdouble left =  -1* diam;
	GLdouble right =  + diam;
	GLdouble bottom=  -1* diam;
	GLdouble top =  diam;
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
    glOrtho(left, right, bottom, top, -100.0, 100.0);
	//glOrtho(left, right, bottom, top, zNear, zFar);
}

void OGLRenderer::drawLines(){
    if(_volPaths.empty())return;
    //glLineWidth(2.0);
    //glBlendFunc(GL_ONE, GL_ONE);
    volumePath *vp;
    for (int i=0; i<_volPaths.size(); i++) {
        vp=_volPaths[i];
        if (vp->points.size()==1) {
            glBegin(GL_POINT);
            glColor4f(1.0,0.,0.,vp->points[0][3]);
            glVertex3f(vp->points[0][0], vp->points[0][1], vp->points[0][2]);
            glEnd();
            continue;
        }
        glBegin(GL_LINE_STRIP);
        for (int j=0; j<vp->points.size(); j++) {
            glColor4f(TREECOLOR,vp->points[j][3]);
            glVertex3f(vp->points[j][0], vp->points[j][1], vp->points[j][2]);
        }
        glEnd();
        glPointSize(3.0);
        glColor4f(1.0f, 1.0f, 0.,0.8);
        glBegin(GL_POINTS);
        for (int j=0; j<vp->points.size(); j++) {
            glVertex3f(vp->points[j][0], vp->points[j][1], vp->points[j][2]);
        }
        glEnd();
    }
}

//Default constructor
OGLRenderer::OGLRenderer(double sceneDiametr){
    defaultScale=1.;
	windowWidth=512;
	windowHeight=512;
	_maxDepthToRender=99999999;

	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition (glutGet(GLUT_SCREEN_WIDTH)/2 - windowWidth/2,
                            glutGet(GLUT_SCREEN_HEIGHT)/2 - windowHeight/2);
	
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA );
	
	glutCreateWindow("Shadow Map Renderer");
	double diam=sceneDiametr;
	double zNear = 0.01;
	double zFar = zNear + diam*5;
	GLdouble left =  -1* diam;
	GLdouble right =  + diam;
	GLdouble bottom=  -1* diam;
	GLdouble top =  diam;
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(left, right, bottom, top, -10000, 10000);
	
	//gluPerspective(0.5 , 1 , zNear ,  zFar );
	
	setupGLUTCallbacks();

}

//my own brute force sleep
void sleep(){
	for (int i=0; i<100000000; i++) {
		//sin(0.34567);
	}
}
void OGLRenderer::setBBoxes(vector<BBox*> bBoxes){
	_bBoxes=bBoxes;
}
//void OGLRenderer::setRoot(Node& root){
//	_root=&root;
//}
void OGLRenderer::drawBox(float min[3],float max[3]){
	float v1[3]={min[0],min[1],min[2]};
	float v2[3]={max[0],min[1],min[2]};
	float v3[3]={min[0],max[1],min[2]};
	float v4[3]={max[0],max[1],min[2]};
	
	float v5[3]={min[0],min[1],max[2]};
	float v6[3]={max[0],min[1],max[2]};
	float v7[3]={min[0],max[1],max[2]};
	float v8[3]={max[0],max[1],max[2]};
	
	//drawing box with lines
	glBegin(GL_LINE_STRIP);
	glVertex3f(v1[0], v1[1], v1[2]);  // V0
	glVertex3f(v2[0], v2[1], v2[2]);  // V1
	glVertex3f(v4[0], v4[1], v4[2]);  // V2
	glVertex3f(v3[0], v3[1], v3[2]);  // V1
	glVertex3f(v1[0], v1[1], v1[2]);  
	
	glVertex3f(v5[0], v5[1], v5[2]); 
	glVertex3f(v6[0], v6[1], v6[2]); 
	glVertex3f(v8[0], v8[1], v8[2]); 
	glVertex3f(v7[0], v7[1], v7[2]); 
	glVertex3f(v5[0], v5[1], v5[2]); 
	glEnd();
	glBegin(GL_LINES);
	glVertex3f(v3[0], v3[1], v3[2]); 
	glVertex3f(v7[0], v7[1], v7[2]);
	
	glVertex3f(v4[0], v4[1], v4[2]);
	glVertex3f(v8[0], v8[1], v8[2]); 
	
	
	glVertex3f(v2[0], v2[1], v2[2]);
	glVertex3f(v6[0], v6[1], v6[2]); 
	glEnd();
}
void OGLRenderer::drawBBoxes(){
	glColor3f(TREECOLOR);
	if (_bBoxes.size()==0) {
		return;
	}
	int ds=_bBoxes[0]->getNumDimensions();
	float min[3];
	float max[3];
	for (int i=0; i<_bBoxes.size(); i++) {
		for (int dim=0; dim<ds; dim++) {
			min[dim]=_bBoxes[i]->getMin(dim);
			max[dim]=_bBoxes[i]->getMax(dim);
		}
		drawBox(min, max);
	}
}

void OGLRenderer::setPoints(vector<float*> points){
    float mnx=1000000000;
    float mny=mnx;
    float mnz=mnx;
    float mxx=-mnx;
    float mxy=-mnx;
    float mxz=-mnx;
    for (int i=0; i<points.size();i++) {
        float x=points[i][0];
        float y=points[i][1];
        float z=points[i][2];
        mnx=(x<mnx)?x:mnx;
        mny=(y<mny)?y:mny;
        mnz=(z<mnz)?z:mnz;
        
        mxx=(x>mxx)?x:mxx;
        mxy=(y>mxy)?y:mxy;
        mxz=(z>mxz)?z:mxz;
    }
    defaultScale=sqrtf(mnx*mnx+mny*mny+mnz*mnz)+sqrtf(mxx*mxx+mxy*mxy+mxz*mxz);
    defaultScale=1./defaultScale;
//    for (int i=0; i<points.size();i++) {
//        points[i][0]*=defaultScale;
//        points[i][1]*=defaultScale;
//        points[i][2]*=defaultScale;
//    }
	_origPoints=points;
}

void OGLRenderer::drawPoints(){

    if (_origPoints.size()==0) {
        return;
    }
	glPointSize(6.0f);
	int sz=min((int)_origPoints.size(), 100000);
	glColor4f(POINTCOLOR,0.3);
	glBegin(GL_POINTS);
	for (int i=0; i<sz; i++) {
		glVertex3f(_origPoints[i][0], _origPoints[i][1], _origPoints[i][2]);
	}
	glEnd();
//	stack< Node* > st;
//
//	st.push(_root);
//	while (st.size()!=0) {
//		Node* nd=st.top();
//		st.pop();
//		vector<float*> dt=nd->getData();
//		glBegin(GL_POINTS);
//		for (int i=0; i<dt.size(); i++) {
//			
//			glVertex3f(dt[i][0], dt[i][1], dt[i][2]);
//		}
//		glEnd();
//		if (nd->isLeaf()) {
//			
//		}else {
//			st.push(nd->getLeft());
//			st.push(nd->getRight());
//		}
//		
//	}
	//traverseTree(*_root);

}
//renders the whole scene
void OGLRenderer::drawScene(void) {
	
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glClear( GL_COLOR_BUFFER_BIT );
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    
    
    
	//glEnable(GL_DEPTH_TEST);
	_depth=0;
	// drawBBoxes();
    drawLines();
	drawPoints();
    
	//debugTree(*_root, *_bBoxes[0]);
	//printf("\nNumber of poits is: %d",_depth);
	glutSwapBuffers();
}

void OGLRenderer::render(){
	glutMainLoop();
}

void OGLRenderer::idleFunc() {
	glutPostRedisplay();
}

void OGLRenderer::handleKeys(unsigned char key, int x, int y) {
	switch (key) {
		case 27:
			//finalizeCUDA();
			exit(0);
		case 'l':
			_maxDepthToRender=0;
			break;
        case 'r':
            this->resetSceneCamera();
		case '+':
            this->translate_z+=55.5;
            this->CameraRefresh();
			_maxDepthToRender+=1;
			break;
		case '-':
            this->translate_z-=55.5;
            this->CameraRefresh();
			if (_maxDepthToRender>=0) {
				_maxDepthToRender-=1;
			}
			break;
		case 'i':
			_maxDepthToRender=99999999;
			break;


	}
}

//----------------------------------------------------------------------
void OGLRenderer::handleMouse(int button, int state, int x, int y)
{
    //handle mouse interaction for rotating/zooming the view
    if (state == GLUT_DOWN) {
        mouse_buttons |= 1<<button;
    } else if (state == GLUT_UP) {
        mouse_buttons = 0;
    }
	
    mouse_old_x = x;
    mouse_old_y = y;
}


//----------------------------------------------------------------------
void OGLRenderer::cameraMotion(int x, int y)
{
    //handle the mouse motion for zooming and rotating the view
    float dx, dy;
    dx = x - mouse_old_x;
    dy = y - mouse_old_y;
	
    if (mouse_buttons & 1) {
        rotate_x += dy * 0.2;
        rotate_y += dx * 0.2;
    } else if (mouse_buttons & 4) {
        translate_z += dy * 0.1;
    }
	
    mouse_old_x = x;
    mouse_old_y = y;
	
    // set view matrix
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    float sc=defaultScale*translate_z;
    glScalef(sc, sc, sc);
    //glTranslatef(0.0, 0.0, translate_z);
    glRotatef(rotate_x, 1.0, 0.0, 0.0);
    glRotatef(rotate_y, 0.0, 1.0, 0.0);
}

void OGLRenderer::CameraRefresh(){
    // set view matrix
   // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    float sc=defaultScale*translate_z;
    glScalef(sc, sc, sc);
    //glTranslatef(0.0, 0.0, translate_z);
    glRotatef(rotate_x, 1.0, 0.0, 0.0);
    glRotatef(rotate_y, 0.0, 1.0, 0.0);

}

//void OGLRenderer::debugTree(Node& currNode, BBox& currBBox){
//	bool curr=false;
//	if (++_depth>=_maxDepthToRender) {
//		return;
//	}
//	if (_depth==_maxDepthToRender-1) {
//		curr=true;
//	}
//	
//	
//	int ds=_bBoxes[0]->getNumDimensions();
//	float min[3];
//	float max[3];
//
//	for (int dim=0; dim<ds; dim++) {
//		min[dim]=currBBox.getMin(dim);
//		max[dim]=currBBox.getMax(dim);
//	}
//	if (curr) {
//		glColor3f(CURRTREECOLOR);
//	}else {
//		glColor3f(VISITEDTREECOLOR);
//	}
//
//	
//	drawBox(min, max);
//	if (curr) {
//		glColor3f(CURRPOINTCOLOR);
//	}else {
//		glColor3f(VISITEDPOINTCOLOR);
//	}
//	vector<float*> dt=currNode.getData();
//	glBegin(GL_POINTS);
//	for (int i=0; i<dt.size(); i++) {
//		
//		glVertex3f(dt[i][0], dt[i][1], dt[i][2]);
//	}
//	glEnd();
//	if (!currNode.isLeaf()) {
//		vector<BBox*>boxes=currBBox.split(currNode.getAxis(),currNode.getSplitPos());
//		debugTree(*currNode.getLeft(),*boxes[0]);
//		debugTree(*currNode.getRight(),*boxes[1]);
//	}
//}
void OGLRenderer::setupGLUTCallbacks(){
	currentInstance=this;
	glutDisplayFunc(drawCallback);
	glutIdleFunc(idleCallback);
	glutMouseFunc(mouseCallback);
    glutMotionFunc(motionCallback);
	glutIdleFunc(idleCallback);
	glutKeyboardFunc(keyCallback);
}