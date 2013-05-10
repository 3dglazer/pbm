/*
 *  PCloud.cpp
 *  KDTree
 *
 *  Created by System Administrator on 5/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *

 */

#include "PCloud.h"

PCloud::PCloud(){
	_seed= time(NULL);
	_numDimensions=3;
	_setDefaultBoundingBox();
}
PCloud::~PCloud(){
}
void PCloud::_setDefaultBoundingBox(){
	_boundingBox.setNumDimensions(_numDimensions); //default setup
	for (int i=0; i<_numDimensions; i++) {
		_boundingBox.setMin(i, -1);
		_boundingBox.setMax(i, 1);
	}
}
void PCloud::randomizeSeed(){
	_seed= time(NULL);
}
void PCloud::setBoundingBox(BBox const &bBox){
	_boundingBox.setNumDimensions(bBox.getNumDimensions());
	for (int i=0; i<bBox.getNumDimensions(); i++) {
		_boundingBox.setMin(i,bBox.getMin(i));
		_boundingBox.setMax(i,bBox.getMax(i));
	}
	
}
void PCloud::setBoundingBox(BBox* const bBox){
	_boundingBox.setNumDimensions(bBox->getNumDimensions());
	for (int i=0; i<bBox->getNumDimensions(); i++) {
		_boundingBox.setMin(i,bBox->getMin(i));
		_boundingBox.setMax(i,bBox->getMax(i));
	}
	
}
void PCloud::setBoundingBox(float* edgeLength,float* boxCenter,int numDimensions ){
	_boundingBox.setNumDimensions(numDimensions);
	_boundingBox.setNumDimensions(numDimensions); //default setup
	for (int i=0; i<numDimensions; i++) {
		_boundingBox.setMin(i, (boxCenter[i]-(edgeLength[i]/2.0)));
		_boundingBox.setMax(i, (boxCenter[i]+(edgeLength[i]/2.0)));
	}
}

void PCloud::generateFlatDistribution(vector<float*>& _points,int numPoints,int numDimensions){
	if (numDimensions!=_numDimensions) {
		_numDimensions=numDimensions;
		_setDefaultBoundingBox();
	}
	_numDimensions=numDimensions;
	srand ( _seed );
	//vector prealocation
	_points.reserve(numPoints);
	float q=0;
	for (int id=0; id<numPoints; id++) {
		float* point=new float[_numDimensions];
		for (int i=0; i<_numDimensions; i++) {
			//bounding box length is important for generating points evenly in each dimension
			q=_boundingBox.getMax(i)-_boundingBox.getMin(i);
			point[i]=float(q)*(float(rand()) / float(RAND_MAX))+_boundingBox.getMin(i);
		}
		_points.push_back(point);
	}
}
void PCloud::setSeed(int newSeed){
	_seed=newSeed;
}
int PCloud::getSeed(){
	return _seed;
}
float  PCloud::_randGauss(){
	float x1, x2, w, y1;// y2;
	
	do {
		x1 = 2.0 * _randF() - 1.0;
		x2 = 2.0 * _randF() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );
	
	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;
	//y2 = x2 * w;
	return y1;
}

float PCloud::_randF(){
	return (float(rand()) / float(RAND_MAX));
}
float PCloud::_randExp(){
	float u =_randF();
	//the minus one is a rate prarametr of the exponential function;
	u=1-u;
	u=log(u);
	u=u/(-1);
	return u;
}
void PCloud::generateGaussianDistribution(vector<float*>& _points,int numPoints,int numDimensions){
	if (numDimensions!=_numDimensions) {
		_numDimensions=numDimensions;
		_setDefaultBoundingBox();
	}

	srand ( _seed );;
	//vector prealocation
	_points.reserve(numPoints);
	float q=0;
	for (int id=0; id<numPoints; id++) {
		float* point=new float(_numDimensions);
		for (int i=0; i<_numDimensions; i++) {
			//bounding box length is important for generating points evenly in each dimension
			q=_boundingBox.getMax(i)-_boundingBox.getMin(i);
			point[i]=(q)*0.1*(_randGauss())+_boundingBox.getMin(i);
		}
		_points.push_back(point);
	}
}

void PCloud::generateScatteredGaussian(vector<float*>& _points,int numPoints, int numNests,int numDimensions){
	if (numDimensions!=_numDimensions) {
		_numDimensions=numDimensions;
		_setDefaultBoundingBox();
	}
	_numDimensions=numDimensions;
	srand ( _seed );
	//vector prealocation
	_points.reserve(numPoints);
	int pointsInNest=numPoints/numNests;
	vector<float*> sampleNest;
	generateGaussianDistribution(sampleNest, pointsInNest, _numDimensions);
	for (int i=0; i<sampleNest.size(); i++) {
		for (int  dim=0; dim<_numDimensions; dim++) {
			sampleNest[i][dim]=(float)(sampleNest[i][dim]);
		}
	}
	vector<float*>nestCenters;
	generateFlatDistribution(nestCenters, numNests, numDimensions);
	for (int n=0; n<numNests; n++) {
		for (int i=0; i<sampleNest.size(); i++) {
			float *point=new float[numDimensions];
			for (int dim=0; dim<numDimensions; dim++) {
				point[dim]=sampleNest[i][dim]*0.5+nestCenters[n][dim]*1.5;
			}
			_points.push_back(point);
		}
	}
}

void PCloud::generateUniformGrid(vector<float*>& _points,int numPoints,int numDimensions){
	if (numDimensions!=_numDimensions) {
		_numDimensions=numDimensions;
		_setDefaultBoundingBox();
	}
	_numDimensions=numDimensions;
	srand ( _seed );
	//vector prealocation
	_points.reserve(numPoints);
	
	
//	float step=(float)pow(( double )volume/(double)numPoints,( double )1/((double)_numDimensions));
//	float* nPInDim=new float[_numDimensions];
	
	
	vector<float*> coordinatesInEachDim;
	int nPInD=ceil(pow((double)numPoints, (double)(1/(double)numDimensions)));
	float* step=_boundingBox.getEdgeLengths();
	for (int i=0; i<numDimensions; i++) {
		step[i]=step[i]/(float)nPInD;
	}
	for(int i =0;i<_numDimensions;i++){
		float* coords=new float[nPInD];
		int id=0;
		for (float j=_boundingBox.getMin(i); j<=_boundingBox.getMax(i); j+=step[i]) {
			coords[id]=j;
			id++;
		}
		coordinatesInEachDim.push_back(coords);
	}
	
	vector<float*> points;
	//pro kazdej krok vynasobit vsechny body ktery uz existujou a pridat do jejich dimenze prislusny cislo;
	//projdu existujici body pridam jim do prislusne dimenze koordinaty 0 prvku. (pokud zadne body neexistujou rovnou vytvorim body jejichz 0 dim je plna koordinat y coordinatesInEachDim)
	//pak pocet zbyvajicich kombinci krat vytvorim kopii existujicich bodu a do aktualni dimenze jim hodim spravnou koordinatu
	//prvni bod v poli
	float* point= new float[_numDimensions];
	for (int i; i<_numDimensions; i++) {
		point[i]=coordinatesInEachDim[i][0];
	}
	points.push_back(point);
	bool quit=false;
	for (int i=0; i<_numDimensions; i++) {
		if (quit) {break;}
		int currPointCount=points.size();
		//nejdriv vsem existujicim nastavim prislusnou koordinatu
		for (int id=0; id<currPointCount; id++) {
			points[id][i]=coordinatesInEachDim[i][0];
		}
		if (quit) {break;}
		//projdu zbytek koordinat v aktualni dimenzi a skopiruju vsechny jiz existujici body a nastavim jim prislusnou dimenzi na spravnou hodnotu 
		for (int p=1; p<nPInD;p++) {
			if (quit) {break;}
			for (int id=0; id<currPointCount; id++) {
				if (points.size()>=numPoints) {quit=true;}
				if (quit) {break;}
				float* pnt=new float[_numDimensions];
				for (int pc=0; pc<_numDimensions; pc++) {
					pnt[pc]=points[id][pc];
				}
				float* pom=coordinatesInEachDim[i];///wtf
				pnt[i]=pom[p];
				points.push_back(pnt);
			}
			if (quit) {break;}
		}
		if (quit) {break;}
	}
	_points=points;
}
