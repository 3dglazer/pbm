/*
 *  KDTree.h
 *  KDTree
 *
 *  Created by System Administrator on 5/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef KDTREE_H
#define KDTREE_H
#include <vector>
#include <deque>
#include <stack>
#include <queue>
#include "BBox.h"
#include "PQNode.h"
#include "Node.h"
#include "NaivePQueue.h"
#include "NaivePointQueue.h"
#include <math.h>
#include "NodeStack.h"
#include "DataDummy.h"
using namespace std;
/// Structure used to compare DataDummy structure returns true if < (used in stl priority queue)
struct DataDummy_Less_compare : binary_function<DataDummy*, DataDummy*, bool> 
{ 
    bool operator()(const DataDummy* d1, const DataDummy* d2) const 
    { 
        return (d1->value < d2->value); 
    } 
}; 

/// Structure used to compare DataDummy structure returns true if > (used in stl priority queue)
struct DataDummy_Greater_compare : binary_function<DataDummy*, DataDummy*, bool> 
{ 
    bool operator()(const DataDummy* d1, const DataDummy* d2) const 
    { 
        return (d1->value > d2->value); 
    } 
}; 


///Serves for comparing pointers to PQNode object instances, without this function priority queue doesnt work
struct PQNode_compare : binary_function<PQNode*, PQNode*, bool> 
{ 
    bool operator()(const PQNode* p1, const PQNode* p2) const 
    { 
        return (p1->getDist() > p2->getDist()); 
    } 
}; 
/**
	\brief KDTree is class which defines the KDTree.
 */
class KDTree {
public:
	/// class responsible for building and performing querries on kdtree.
	KDTree();
	~KDTree();
	///Retruns true if the kdtree doesnot contain any nodes.
	bool isEmpty();
	///Returns root of the KDTree
	Node* getRoot();
	///Returns value of the traversation profiling variable. Used to measure how many kdtree nodes were visited.
	int getDepthProfiler();
	///Returns all nodes in the form of bounding boxes, this is used for visualization and debuging purpouses.
	vector<BBox*> getBBoxes();
	///Returns bounding box of all the points. (Bounding box of the root node.)
	BBox* getBBox();
	///Builds the kdtree.
	void build(bool recursion);
	/**
	 \brief Method used to perform rectangular range searches.
	 \param range hyper rectangle defining search area.
	 \param recursion if set true recursive method is used else nonrecursive stack based.
	  \return vector containing all points in the range.
	 */
	vector<float*> rangeSearch(BBox const &range,bool recursion);
	/**
	 \brief Method used to perform circular range search.
	 \param center float array containing center coordinates of the hyper circle.
	 \param radius is a radius of the hyper cirle.
	 \param recursion if set true recursive method is used else nonrecursive priority queue version.
	 \return vector containing all points in the range.
	 */
	vector<float*> circularSearch(float* center,float radius,bool recursion);
	/**
	 \brief Method used to perform k-nearest neighbour search querries.
	 \param q float array containg coordinates of the query point.
	 \param k count of nearest neighbours wanted.
	 \param recursion if set true recursive method is used else nonrecursive priority queue version.
	 \return vector containing k nearest neighbours to query point q.
	 */
	vector<float*> nnSearch(float* q, int k,bool recursion);
	/**
	 \brief This function sets the points. 
	 */
	void setData(vector<float*>& data,int dimensionality );
	/// Sets the division criteria, how many points should be at maximum at leaf nodes.
	void setMaxPointsInLeaf(int max);
	/// Sets the type of priority queue 1 is naive priority queue and 0 is stl_vector priority queue.
	void setPriorityQueueType(int pqType);
	void setStackType(int stType);
	/// Constant specifing priority queue type
	static const int STL_VECTOR_PQUEUE = 0;
	/// Constant specifing priority queue type
	static const int NAIVE_PQUEUE= 1;
	static const int  STL_STACK=0;
	static const int  NODE_STACK=1;
private:
	/// rangeSearch using stack.
	void _rangeSearch(vector<float*>& data,BBox const &range,Node const &rt);
	///	range search implemented as recursion.
	void _recursiveRangeSearch(vector<float*>& data,BBox const &range,Node const &rt);
	///Nonrecursive implmenentation of nearest neighbour search using priority queue. The type of pq can be set with setPriorityQueueType function.
	void _nnSearch(); //nonrecursive implementation with priority queue
	
	/**
	 \brief Nonrecursive implmenentation of circular search using priority queue. The type of pq can be set with setPriorityQueueType function.
	 \param data this vector will be filled with points contained in the specified circular region.
	 */
	void _circularSearch(vector<float*>& data); //nonrecursive implementation with priority queue
	/**
	 \brief recursive implementation of nearest neighbour search.
	 \param u is current node to be examined.
	 \param rd is distance to childs of currently examined node.
	 */
	void _recursiveSearchNN(Node const &u,float rd );
	/**
	 \brief recursive implementation of circular range search.
	 \param data this vector will be filled with points contained in the specified circular region.
	 \param u is current node to be examined.
	 \param rd is distance to childs of currently examined node.
	 */
	void _recursiveCircularSearch(vector<float*>& data,Node const &u,float rd ); 
	
	/**
	 \brief nonrecursive implementation of k nearest neighbour search using two priority queues, one for nodes and one for points.
	 */
	void _kNNSearch(vector<float*>&data,int k);
	/**
	 \brief Computes bounging box of the data.
	 */
	void _computeBBox();
	/**
	\brief Recursive function used to build the kdtree.
	 \param data2Split indices of points that should be split.
	 \param axs is an index of axis used to split last node.
	 \param depth contains an information in which level in the kdtree is the current node. This info can be used for choosing the split axis.
	 \param currBBox bounding box of the current node.
	 */
	Node* _build(vector<int> const &data2Split,int axs,int depth, BBox& currBBox);
	///This function is used in recursive kdtree build. It chooses split dimension to split, and fills in dataLeft inices and dataRight indices.
	float _split(vector<int> const &data2Split,vector<int>& dataLeft,vector<int>& dataRight,int& axis, BBox& currBBox);
	///This function counts number of ellements on the left of the split plane and on the right of the split plane. 
	void  _split(int* const data2Split,int &data2splitSize,int& dataLeftSize,int& dataRightSize,int& axis, BBox& currBBox);
	///This function fills dataLeft and dataRight arrays with 
	float _split(int* const data2Split,int &data2splitSize,int* dataLeft,int* dataRight,int &axis,float &splitPos);
	/// Computes Euclidian distance between two multidimensional points.
	inline float _dist2p(const float* a, const float* b);
	/// Computes squared Euclidian distance between two multidimensional points.
	inline float _dist2ps(const float* a, const float* b);
	///Query point for nnSearch also center for circular search
	float* _q; 
	/// Array of offsets for nnSearch.
	float* _off; 
	/// Neares neighbour point.
	float* nnPoint;
	/// Neares neighbour distance and also serves as radius in circular search.
	float nn_dist;
	
	///Variable _depth is used in profiling, measures number of kdtree's nodes visited during searches.
	int _depth;
	/// Root node of the kdtree.
	Node* _root;
	///local copy of the poits seted by setData funcion.
	vector<float*> _data;
	///dimensionality of the kdtree
	int _dimensionality;
	///bounding box of the data;
	BBox* _worldBoundingBox;
	///for debuginng purpouses, stores hyper rectangle for each node
	vector<BBox*> _kdTree;
	///this is kdtree build termination criteria, the maximal point count in the leaf node.
	int _maxDataCount;
	///this is kdtree build termination criteria, the maximal depth of the kdtree.
	int _maxDepth;
	///This sets priority queue type used in searches.
	int _priorityQueType;
	///Sets which stack should be used.
	int _stackType;
	/// My implementation of the stack
	NodeStack nodeSt;
	/// stl implementation of the stack
	stack< Node* > st;
	
	
};
#endif