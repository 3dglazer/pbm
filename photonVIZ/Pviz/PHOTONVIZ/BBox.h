/*
 *  BBox.h
 *  KDTree
 *
 *  Created by System Administrator on 5/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef BBOX_H
#define BBOX_H	
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;
///An implementation of hyper rectangle, you can create any n-dimensional hyper rectagle. Usefull for range querries.
class BBox {
public:
	///Sets minimum value of hyper rectangle in specified dimension.
	void setMin(int dimension, float value);
	///Sets maximum value of hyper rectangle in specified dimension.
	void setMax(int dimension, float value);
	///Returns minimum value of hyper rectangle in specified dimension.
	float getMin(int dimension) const;
	///Returns maximum value of hyper rectangle in specified dimension.
	float getMax(int dimension) const;
	///Retruns float array of lengths off all the edges in all the dimensions.
	float* getEdgeLengths() const;
	///Sets number of dimensions hyper rectangle will hold.
	void setNumDimensions(int dims); // sets number of dimensions
	///Retruns number of dimensions of the hyper rectangle.
	int getNumDimensions() const;
	///Prints all the minimum and maximum values of the hyper rectangle
	void toString() const;
	///Splits the hyper rectangle in the specified dimension in a specified position (worldValue) and returns vector containing two pointers to newly created BBoxes.
	vector<BBox*> split(int dimension, float worldValue) const; //returns two new bounding Boxes
	///Returns int index of the longest edge of the hyper rectangle.
	int getLongest() const;//returns longest dimension of the bounding box
	///Default constructor
	BBox();
	///Constructor with array initialization and sets number of dimensions.
	BBox(int numOfDimensions);
	~BBox();
	
private:
	///World coordinates of the vertices which define the hyper bounding box
	vector<float> _bBox;
	///Number of dimensions
	int _nDims;
};
#endif