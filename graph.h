#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

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
}//end binary heap class
;

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
}
;

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
	BinaryMinHeaps(int size, int thread); //constructor
	~BinaryMinHeaps(); //destructor
	void free(int thread);//theaded destructor
	void siftUp(int nodeIndex, int thread);
	void insert(Nodep* Node, int thread);
	void siftDown(int nodeIndex, int thread);
	Nodep* extractMin(int thread);
	Nodep* extractMintwee(int thread);
}//end binary heaps class
;

class Graph
{
	int maxnodes;
	int maxarcs;
	vector<Link> l;
	vector<Node> n;
	vector<Node> np;
	Graph(int maxnodes,int maxlinks, int maxtimeslots) {};
};

