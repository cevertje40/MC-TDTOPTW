#pragma once
#include "instance.h"

class Sol
{
public:
	vector<vector<Ins::Vertex*>> solution;//solution multiple paths containing a sequence of vertex pointers
	vector<vector<int>> traveltime;// multiple paths containing departuretime-t[pathindex].EDT at each vertex
	vector<vector<int>> max_shift;//local evaluation metric, maximum amount of time each vertex can be shifted forward in time
	vector<vector<int>> action;//0 visit 1 visit and break
	boost::dynamic_bitset<> available;// bitset that states for every vertex if it is still available for inclusion
	vector <int> scores;//score per tour
	vector <int> weight;//weight per tour
	vector <int> volume;//volume per tour
	int score;//sum of path scores
	vector <int> breakindex;//position of break per tour
	//methods
	Sol()
	{

	}
	~Sol()
	{

	}
	Sol(int& maxvertices, int& maxtour, Ins::Vertex* startdepot)
	{
		solution.resize(maxtour);
		traveltime.resize(maxtour);
		max_shift.resize(maxtour);
		scores.resize(maxtour);
		weight.resize(maxtour);
		volume.resize(maxtour);
		action.resize(maxtour);
		breakindex.resize(maxtour);
		for (int t = 0; t < maxtour; ++t)
		{
			solution[t].reserve(maxvertices);//reserves memory
			solution[t].push_back(startdepot);//insert start depot
			traveltime[t].reserve(maxvertices);//reserve memory
			traveltime[t].push_back(0);//insert first travel time
			action[t].reserve(maxvertices);
			action[t].push_back(0);
			max_shift[t].reserve(maxvertices);
			max_shift[t].push_back(0);
			scores[t] = 0;
			weight[t] = 0;
			volume[t] = 0;
			breakindex[t] = -1;
		}
		score = 0;
		available = boost::dynamic_bitset<>(maxvertices);
		available.set();//sets all bits to true
		available[startdepot->index] = false;
	}
};