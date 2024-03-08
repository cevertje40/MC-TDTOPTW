#pragma once
#include "graph.h"

using namespace std;

class Ins//problem instance class that stores all required information f
{
public:
	class Connec
	{
	public:
		//variables
		int from;//vertex index of departing vertex
		int to;// vertex index of arrival vertex
		int dist;//length in km
		int determin;// deterministic time-independent travel time
		vector<double> mu;//for deterministic time-dependent travel time for every timeslot
		vector<double> nu;//for deterministic time-dependent travel time for every timeslot
		//methods
		Connec() {}
		~Connec() {}
		Connec(int& from, int& to, int& dist, int& determin) : from(from), to(to), dist(dist), determin(determin) {}
	};

	class Vertex
	{
	public:
		int id;// order id
		int index;//index in the vertex vector
		int score;//score
		double serv;//service time
		vector<double> LTW;//lower time window for per tour
		vector<double> UTW;//upper time window for per tour
		double weight;//weight
		double volume;//volume
		vector<vector<Vertex*>> nb; //pointer set of neighhours for each day
		vector<boost::dynamic_bitset<>> nbi;// bitset of neighbours for each day
		vector<Connec*> con;// pointer set of connections leaving from the vertex under consideration
		//methods
		Vertex() {}
		~Vertex() {}
		Vertex(int id, int index, int score, vector<double> LTW, double serv, vector<double> UTW, double weight, double volume) : id(id), index(index), score(score),serv(serv), LTW(LTW), UTW(UTW), weight(weight), volume(volume) {}
	};

	class Tour
	{
	public:
		int id;
		int index;//dayindex, to obtain corresponding time windows
		Vertex* startv;//vertex pointer to start depot
		Vertex* endv;//vertex pointer to end depot
		int EDT;//earliest departure time
		int LAT;//latest arrival time
		int T_max;//maximum allowable travel time
		int W_max;//maximum allowable weight
		int V_max;//maximum allowable volume
		Tour() {}
		~Tour() {}
		Tour(int id, int index, Vertex* startv, Vertex* endv, int EDT, int LAT, int T_max, int W_max, int V_max) : id(id), index(index), startv(startv), endv(endv), EDT(EDT), LAT(LAT), T_max(T_max), W_max(W_max), V_max(V_max) {}
	};

	int maxvertices;
	int maxtours;
	double breakstart;
	double breakend;
	double breakdur;
	vector<Vertex> v;//vertex objects
	vector<Connec> c;//connection objects
	vector<Tour> t;//tour objects

	Ins(string filename);//construct instance by reading file
	void construct_time_independent_traveltime(Graph& graph);
	void construct_time_dependent_traveltime(Graph& graph);

};
