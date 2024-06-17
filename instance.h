#pragma once
#include "graph.h"

using namespace std;

class Ins//problem instance class that stores all required information f
{
public:
	struct MCTDTOPTW 
	{
		string path;
		string name;
		
	};

	struct CTOP {
		string path;
		string name;
	};

	class Connec
	{
	public:
		//variables
		int from;//vertex index of departing vertex
		int to;// vertex index of arrival vertex
		double determin;// deterministic time-independent travel time
		vector<double> mu;//for deterministic time-dependent travel time for every timeslot
		vector<double> nu;//for deterministic time-dependent travel time for every timeslot
		//methods
		Connec() {}
		~Connec() {}
		Connec(int& from, int& to, double& determin) : from(from), to(to), determin(determin) {}
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
		double EDT;//earliest departure time
		double LAT;//latest arrival time
		double T_max;//maximum allowable travel time
		double W_max;//maximum allowable weight
		double V_max;//maximum allowable volume
		Tour() {}
		~Tour() {}
		Tour(int id, int index, Vertex* startv, Vertex* endv, double EDT, double LAT, double T_max, double W_max, double V_max) : id(id), index(index), startv(startv), endv(endv), EDT(EDT), LAT(LAT), T_max(T_max), W_max(W_max), V_max(V_max) {}
	};

	int maxvertices;
	int maxtours;
	double breakstart;
	double breakend;
	double breakdur;
	double maxscore;//total score over all vertices
	vector<Vertex> v;//vertex objects
	vector<Connec> c;//connection objects
	vector<Tour> t;//tour objects
	default_random_engine engine;

	
	Ins(struct MCTDTOPTW);//construct instance by reading file
	Ins(struct CTOP);//convert CTOP instance

	void construct_time_independent_traveltime(Graph& graph);
	void construct_time_dependent_traveltime(Graph& graph);
	void read_time_independent_traveltime();
	void read_time_dependent_traveltime();
	void create_neighbourhood(string path, string name,int amnt_nb);
	void read_neighbourhood(string path,string name);
	
	/**acces of c object methods*/
	int find_t(double time);
	double travel_time(Connec* c, double start);
	double arrival_time(Connec* c, double start);
	double departure_time(Connec* c, double arrivaltime);
};
