/*
 *  PCloud.h
 *  KDTree
 *
 *  Created by System Administrator on 5/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  This is a point cloud generator. Can generate Flat distribution (using standart rand()),
 *	uniform grid (n^1/3 is not a integer, one side will not be compleated.),
 *	gaussian distribution, and scattered gaussian distribution.
 */

#ifndef PCLOUD_H	
#define PCLOUD_H
#include "BBox.h"
#include <vector>
#include <stdlib.h>
#include <math.h>
using namespace std;
/// Class which can generate point clouds with variouse distributions.
class PCloud {
public:
	PCloud();
	~PCloud();
	///Sets approximate bounding box that the points will try to obey.
	void setBoundingBox(BBox const &bBox);
	///Sets approximate bounding box that the points will try to obey.
	void setBoundingBox(BBox* const bBox);
	/**
	 \brief Sets approximate bounding box.
	 \param edgeLength float array containing edge lengths in all dimensions.
	 \param boxCenter float array containing position of the center of the box.
	 \param numDimensions sets number of dimensions of the box.
	 */
	void setBoundingBox(float* edgeLength,float* boxCenter,int numDimensions );
	/**
	 \brief Generates points with gaussian distribution.
	 \param _points Vector which will be filled with generated points (n-dimensional float arrays).
	 \param numPoints Number of poits to generate.
	 \param numDimensions Sets dimensionality of the generated points.
	 */
	void generateGaussianDistribution(vector<float*>& _points,int numPoints,int numDimensions);
	
	/**
	\brief  Generates points with flat distribution, using standart rand() function.
	\param _points Vector which will be filled with generated points (n-dimensional float arrays).
	\param numPoints Number of poits to generate.
	\param numDimensions Sets dimensionality of the generated points.
	 */
	void generateFlatDistribution(vector<float*>& _points,int numPoints,int numDimensions);
	/**
	 \brief Generates points with scattered gaussian distribution.
	 \param _points Vector which will be filled with generated points (n-dimensional float arrays).
	 \param numPoints Number of poits to generate.
	 \param numNests Number of point nests (clusters of points with gaussian distribution).
	 \param numDimensions Sets dimensionality of the generated points.
	 */
	void generateScatteredGaussian(vector<float*>& _points,int numPoints, int numNests,int numDimensions);
	/**
	 \brief Generates uniform grid which will contains numPoints number of points with dimensionality specified in numDimensions parametr.
	 */
	void generateUniformGrid(vector<float*>& _points,int numPoints,int numDimensions);
	/// Sets seed for random generator.
	void setSeed(int newSeed);
	/// Sets new random seed, based on the current computer time.
	void randomizeSeed();
	/// Returns seed used for generating random numbers.
	int	 getSeed();
private:
	/**
		\brief Returns random float number between 0-1 using random generator witch gaussian distribution.
	 */
	float _randGauss();
	/**
	  \brief Returns random float number between 0-1 using random generator witch flat distribution.
	 */
	float _randF();
	///Returns random number with gaussian (exponential) distribution.
	float _randExp();
	///Creates deffault bounding box, with corners at (-1,1) in all dimensions.
	void _setDefaultBoundingBox();
	///Number of dimensions of the generated points.
	int _numDimensions;
	///Approximate bounding box of the points generated, will not be strictly obeyed!
	BBox _boundingBox;
	///Seed for random number generator.
	int _seed;
};

#endif