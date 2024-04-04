#pragma once
#include "instance.h"

class Sol
{
public:
	vector<vector<Ins::Vertex*>> solution;//solution multiple paths containing a sequence of vertex pointers
	vector<vector<double>> traveltime;// multiple paths containing departuretime-t[pathindex].EDT at each vertex
	vector<vector<double>> max_shift;//local evaluation metric, maximum amount of time each vertex can be shifted forward in time
	vector<vector<int>> action;//0 visit, 1 visit and break
	boost::dynamic_bitset<> available;// bitset that states for every vertex if it is still available for inclusion
	vector <int> scores;//score per tour
	vector <double> weight;//weight per tour
	vector <double> volume;//volume per tour
	int score;//sum of path scores
	vector <int> breakindex;//position of break per tour
	Ins* ins;//pointer to instance object
	//methods
	Sol(){}
	~Sol(){}
	Sol(Ins& ins);
	void reset();
	void check();
	void update_traveltime(int tour, int start, int end);
	void calc_maxshift(int tour, int start, int end);
	void calc_maxshift();
	void insertvertex(int tour, Ins::Vertex* insertion, int position);
	void removevertex(int tour, int position);
};