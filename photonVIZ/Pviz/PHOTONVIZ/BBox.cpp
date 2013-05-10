/*
 *  BBox.cpp
 *  KDTree
 *
 *  Created by System Administrator on 5/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights resized.
 *
 */

#include "BBox.h"
BBox::BBox(){
	_nDims=0;
}

float* BBox::getEdgeLengths() const{
	float* ret=new float[_nDims];
	for (int i=0; i<_nDims; i++) {
		ret[i]=fabs(float(_bBox[(2-1)+1]-_bBox[(2-1)+2]));
	}
	return ret;
}
BBox::BBox(int numOfDimensions){
	_nDims=numOfDimensions;
	_bBox.resize(_nDims*2);
}
BBox::~BBox(){
}
int BBox::getNumDimensions() const{
	return _nDims;
}
void BBox::setNumDimensions(int dims){
	_nDims=dims;
	_bBox.resize(_nDims*2);
}
void BBox::setMin(int dimension,float value){
	if (dimension>=_nDims) {
		//printf("\nBounding Box Error:: dimension %d out of bounds\n",dimension);
	}else {
		_bBox[(dimension*2-1)+1]=value;
	}
}
void BBox::setMax(int dimension,float value){
	if (dimension>=_nDims) {
		//printf("\nBounding Box Error:: dimension %d  out of bounds\n",dimension);
	}else {
		_bBox[(dimension*2-1)+2]=value;
	}
}
float BBox::getMin(int dimension) const{
	if (dimension>=_nDims) {
		//printf("Bounding Box Error:: dimension out of bounds");
		return 0;
	}else {
		return _bBox[(dimension*2-1)+1];
	}
}
float BBox::getMax(int dimension) const{
	if (dimension>=_nDims) {
		//printf("Bounding Box Error:: dimension out of bounds");
		return 0;
	}else {
		return _bBox[(dimension*2-1)+2];
	}
}
int BBox::getLongest() const{
	float max=fabs(float(_bBox[0]-_bBox[1]));
	int maxId=0;
	float tmp=0;
	for (int i=1; i<_nDims; i++) {
		tmp=fabs(float(_bBox[(i*2-1)+1]-_bBox[(i*2-1)+2]));
		if (tmp>max) {
			max=tmp;
			maxId=i;
		}
	}
	return maxId;
}
 void BBox::toString() const{
	printf("||BBox:");
	for (int i=0; i<_nDims; i++) {
		printf("[%4.2f,%4.2f]",getMin(i),getMax(i));
	}
	printf("BBox||");
}
vector<BBox*> BBox::split(int dimension, float worldValue) const{
	if (dimension>_nDims && _bBox.size()==0) {
		vector<BBox*> b;
		return b;
	}else {
		vector<BBox*> bBoxes;
		bBoxes.push_back(new BBox(_nDims));
		bBoxes.push_back(new BBox(_nDims));
		//memcpy(left,this,sizeof(BBox));
		for (int i=0; i<_nDims; i++) {
			bBoxes[0]->setMin(i, getMin(i));
			bBoxes[0]->setMax(i, getMax(i));
			bBoxes[1]->setMin(i, getMin(i));
			bBoxes[1]->setMax(i, getMax(i));
		}
		bBoxes[0]->setMax(dimension, worldValue);
		//memcpy(right,this,sizeof(BBox));
		bBoxes[1]->setMin(dimension, worldValue);
		return bBoxes;
	}
}