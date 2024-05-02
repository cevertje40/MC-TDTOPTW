#pragma once
#include "instance.h"



class Sol
{
public:
	class Tour
	{
		public:
		int index;
		vector <Ins::Vertex*> seq;//sequence of vertex pointers
		vector<double> deptime;//departuretime-EDT at each vertex
		vector<double> max_shift;//local evaluation metric, maximum amount of time each vertex can be shifted forward in time
		vector<int> action;//0 visit, 1 break and visit
		int score;//total score of tour
		double weight;//weight per tour
		double volume;//volume per tour
		int breakindex;//position of break in tour
	};
	vector<Tour>tours;//solution consist of collection of tours
	vector<int>tourindex;//random tour index
	boost::dynamic_bitset<> available;// bitset that states for every vertex if it is still available for inclusion
	int score;//sum of all tour scores
	Ins* ins;//pointer to instance object
	//methods
	Sol(){}
	~Sol(){}
	Sol(Ins& ins);
	void reset();
	void check();
	void update_traveltime(int tour, int start, int end);
	void update_traveltime_break(int tour, int start, int end);
	void update_maxshift(int tour, int start, int end, double arrivaltime);//partial update within a tour
	void calc_maxshift(Sol::Tour &tour);//for specific tour of solution
	void insertvertex(Sol::Tour &tour, Ins::Vertex* candidate, int position);
	void replacevertex(Sol::Tour &tour, Ins::Vertex* candidate, int position, bool updatebreak=false);
	void removevertex(Sol::Tour &tour, int position);
	void swapvertex(Sol::Tour &tour, int i, int j);//assumption i < j
	void optvertices(Sol::Tour &tour, int i, int j);//assumption i < j
	friend ostream& operator<<(ostream& output, Sol& sol);
	
};