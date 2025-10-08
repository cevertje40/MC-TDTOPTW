#pragma once
#include "graph.h"

class Ins//problem instance class that stores all required information
{
public:
	struct MCTDTOPTW 
	{
		std::string path;
		std::string name;
		
	};

	struct CTOP {
		std::string path;
		std::string name;
	};

	class Connec
	{
	public:
		//variables
		int from;//vertex index of departing vertex
		int to;// vertex index of arrival vertex
		double determin;// deterministic time-independent travel time
		std::vector<double> mu;//for deterministic time-dependent travel time for every timeslot
		std::vector<double> oneplusmu;//1+mu for faster computation
		std::vector<double> nu;//for deterministic time-dependent travel time for every timeslot
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
		std::vector<double> LTW;//lower time window for per tour
		std::vector<double> UTW;//upper time window for per tour
		double weight;//weight
		double volume;//volume
		std::vector<std::vector<Vertex*>> nb; //pointer set of neighhours for each day
		std::vector<boost::dynamic_bitset<>> nbi;// bitset of neighbours for each day
		std::vector<Connec*> con;// pointer set of connections leaving from the vertex under consideration
		//methods
		Vertex() {}
		~Vertex() {}
		Vertex(int id, int index, int score, std::vector<double> LTW, double serv, std::vector<double> UTW, double weight, double volume) : id(id), index(index), score(score),serv(serv), LTW(LTW), UTW(UTW), weight(weight), volume(volume) {}
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
	std::vector<Vertex> v;//vertex objects
	std::vector<Connec> c;//connection objects
	std::vector<Tour> t;//tour objects

	
	Ins(struct MCTDTOPTW);//construct instance by reading file
	Ins(struct CTOP);//convert CTOP instance

	void construct_time_independent_traveltime(Graph& graph);
	void construct_time_dependent_traveltime(Graph& graph);
	void read_time_independent_traveltime();
	void read_time_dependent_traveltime();
	void create_neighbourhood(std::string path, std::string name);
	void read_neighbourhood(std::string path, std::string name);

	void alter_instance();
	void unalter_instance();
	
	/**acces of c object methods*/
	inline int find_t(double time)
	{
		int t = (int)((time - time_periods[0]) / 0.25);//when you change the time unit this has to change too
		return std::min(55, t);
	}

	inline double travel_time(Connec* c, double start)
	{
		int t = find_t(start);
		return std::fma(c->nu[t] - 1.0, start, c->mu[t]);
		/*
		int t = find_t(start);
		double traveltime = c->nu[t] + start * c->mu[t];
		return traveltime;
		*/
	}

	inline double arrival_time(Connec* c, double start)
	{
		const int t = find_t(start);
		// arr = a + b*start = nu + (1+mu)*start
		return std::fma(c->oneplusmu[t], start, c->nu[t]);
		/*
		int t = find_t(start);
		double arrivaltime = c->nu[t] + (start)*c->mu[t] + start;
		return arrivaltime;
		*/
	}

	inline double departure_time(Connec* c, double arrivaltime)
	{
		int t = find_t(arrivaltime);             // initial guess, not guaranteed correct
		for (;;) {
			// require FIFO: b > 0
			double start = (arrivaltime - c->nu[t]) / c->oneplusmu[t];
			if (start < time_periods[t]) {
				if (t == 0) return time_periods[0];
				--t; continue;
			}
			/*
			if (start >= time_periods[t + 1]) 
			{
				if (t >= 55 - 2) return std::nextafter(time_periods[55 - 1], -INFINITY);
				++t; continue;
			}
			*/
			return start; // landed in the correct bucket
		}
		/*
		int t = find_t(arrivaltime);
		double departuretime = (arrivaltime - c->nu[t]) / (c->oneplusmu[t]);
		while ((time_periods[t] > departuretime) || (departuretime > time_periods[t + 1]))
		{
			--t;
			departuretime = (arrivaltime - c->nu[t]) / (c->oneplusmu[t]);
		}
		return departuretime;
		*/
	}

};
