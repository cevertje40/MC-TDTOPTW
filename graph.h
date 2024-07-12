#pragma once
#include "ConsoleColor.h"
#include <fstream>//input
#include <string.h>
#include <string>
#include <time.h>// cpu time &date
#include <map>
#include <list>
#include "math.h"
#include <algorithm> 
#include <tuple>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
#include <omp.h>
#include <boost/dynamic_bitset.hpp>
#include "time.h"
#include <boost/heap/priority_queue.hpp>
#include <boost/circular_buffer.hpp>


using namespace std;

class Link
{
private:

public:
	int link_id;//id is index+1
	int from;// origin of arc index niet id
	int to;	// TO: destination of arc index niet  id
	double optimaltt;	//in minute
	Link* next;	// LL control
	Link* nextbackward;// LL control
	vector<double> traveltime;//in minute
	vector<double> k;
	vector<double> th;

	Link() : link_id(0), from(0), to(0), optimaltt(0), next(0), nextbackward(0)
	{

	}

	~Link()
	{
	}
};

class Node
{
public:
	int id;//id is index+1
	double longitude;
	double latitude;
	Link* first;
	Link* firstbackward;
	//dijkstra related
	double score;//real score in dijkstra and estimate in a-star
	bool istarget;
	int hdex; // Location(index) in heap
	//function
	Node(const double& KEY = DBL_MAX) : score(KEY), first(0), firstbackward(0), id(0), longitude(0.0), latitude(0.0), istarget(false), hdex(0) {}
	~Node()//destructor
	{
	}
	void addarc(Link* ARC);
	void addprevious(Link* ARC);
};

class BinaryMinHeap
{//binary min heap used in the single threaded dijkstra algorithm
private:
	Node** Nodes;
	int heapSize;
	inline  int getLeftChildIndex(int nodeIndex);
	inline   int getRightChildIndex(int nodeIndex);
	inline   int getParentIndex(int nodeIndex);
public:
	BinaryMinHeap(int size);
	~BinaryMinHeap();
	void siftUp(int nodeIndex);
	void insert(Node* Node);
	void siftDown(int nodeIndex);
	Node* extractMin();//end extractmin
	Node* extractMintwee();//end extractmintwee
};//end binary heap class

class Nodep
{
public:
	int id;
	double longitude;
	double latitude;
	Link* first;
	Link* firstbackward;
	//dijkstra related
	double score[12];//real score in dijkstra estimate in a-star
	int hdex[12]; // Location(index) in heap
	bool istarget[12];
	//function
	Nodep() :first(0), firstbackward(0), id(0), longitude(0.0), latitude(0.0), score{ 0 }, hdex{ 0 }, istarget{ 0 } {}
	~Nodep()//destructor
	{
		//delete first;
		//first=NULL;
		//delete firstbackward;
		//firstbackward=NULL;
	}
	void addarc(Link* ARC);
	void addprevious(Link* ARC);
};

class BinaryMinHeaps
{//binary min heap for parallel computing using in the threaded dijkstra algorithm
private:
	Nodep** Nodes[12];//assuming 12 threads
	int heapSize[12];
	inline  int getLeftChildIndex(int nodeIndex);
	inline   int getRightChildIndex(int nodeIndex);
	inline   int getParentIndex(int nodeIndex);
public:
	BinaryMinHeaps(); // default constructor
	BinaryMinHeaps(int size, int thread); //threaded constructor
	~BinaryMinHeaps(); //destructor
	void free(int thread);//threaded destructor
	void siftUp(int nodeIndex, int thread);
	void insert(Nodep* Node, int thread);
	void siftDown(int nodeIndex, int thread);
	Nodep* extractMin(int thread);
	Nodep* extractMintwee(int thread);
};//end binary heaps class

class Graph
{
private:
	double calculate_mean_d(vector<double>& input);
	double calculate_stdv_pop_d(double mean, vector<double>& input);
	double calculate_stdv_sample(double mean, vector<double>& input);
public:
	int maxnodes;
	int maxlinks;
	
	vector<Link> l;
	vector<Node> n;
	vector<Nodep> np;
	Graph(int maxnodes,int maxlinks);
	double dijkstra_independent(int source, int target);//TI 1 to 1
	vector<double> dijkstra_independent_to_all_threaded(int source, vector<int> targets, int thread);//TI 1 to all thread safe
	double dijkstra_dependent(int source, int target, double currenttime);//TD 1 to 1
	vector<double> dijkstra_dependent_to_all_threaded(int source, vector<int>targets, double currenttime, int thread); //TD 1 to all thread safe
};

