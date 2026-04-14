#pragma once
#include "instance.h"


struct RemovedVertexStats
{
	int removed_count = 0;
	int removed_score_sum = 0;
	double removed_tw_width_sum = 0.0;
	double removed_service_sum = 0.0;
	double removed_weight_sum = 0.0;
	double removed_volume_sum = 0.0;
	double removed_depot_tt_sum = 0.0;
	double removed_position_sum = 0.0; // normalized position in route
};

struct RoutePerfStats
{
	int served_customers = 0;
	int collected_score = 0;

	double route_duration = 0.0;
	double route_time_budget = 0.0;
	double route_time_utilization = 0.0;

	double weight_used = 0.0;
	double weight_capacity = 0.0;
	double weight_utilization = 0.0;

	double volume_used = 0.0;
	double volume_capacity = 0.0;
	double volume_utilization = 0.0;

	double total_waiting_time = 0.0;
	double avg_waiting_time_per_customer = 0.0;

	bool break_taken = false;
	double break_start_time = -1.0;      // elapsed from EDT
	double break_position_norm = -1.0;   // normalized customer position
	bool break_at_end_depot = false;

	RemovedVertexStats removed;
};

struct SolutionPerfStats
{
	int total_score = 0;
	int total_served_customers = 0;
	int total_removed_customers = 0;
	double avg_score_per_served_customer = 0.0;
	double avg_route_duration = 0.0;
	double avg_route_time_utilization = 0.0;
	double avg_weight_utilization = 0.0;
	double avg_volume_utilization = 0.0;
	double avg_waiting_time_per_customer = 0.0;
	double avg_break_start_time = 0.0;
	double avg_break_position_norm = 0.0;
	double pct_break_at_end_depot = 0.0;
	double avg_removed_score = 0.0;
	double avg_removed_tw_width = 0.0;
	double avg_removed_service = 0.0;
	double avg_removed_weight = 0.0;
	double avg_removed_volume = 0.0;
	double avg_removed_depot_tt = 0.0;
	double avg_removed_position = 0.0;
	std::vector<RoutePerfStats> route_stats;
};


class Tour
{
	public:
	Ins* ins;//pointer to instance object
	int index;
	std::vector <Ins::Vertex*> seq;//sequence of vertex pointers
	std::vector<double> deptime;//departuretime-EDT at each vertex
	std::vector<double> max_shift;//local evaluation metric, maximum amount of time each vertex can be shifted forward in time
	std::vector<int> action;//0 visit, 1 break and visit
	int score;//total score of tour
	double weight;//weight per tour
	double volume;//volume per tour
	int breakindex;//position of break in tour
	void update(int start, int end);//keep break fixed and update time and maxshift for positions after start
	void update_break(int newbreakindex);//update solution and position break at input breakindex
	void update_break();//complete update of the tour with optimized break repositioning, only works for feasible tours
	void update_maxshift(int start, int end, double arrivaltime);//update maxshift for positions before end
	void calc_maxshift();//for specific tour of solution
	void insert_vertex(Ins::Vertex* candidate,int position);
	void remove_vertex(int position);
	void remove_vertices(int position1, int position2);
	void replace_vertex(Ins::Vertex* candidate, int position);
	void replace_vertex(Ins::Vertex* candidate, int position,int breakindex);
	void opt_vertices(int i, int j);//assumption i < j
	void swap_vertices(int i, int j);//assumption i < j
	bool check();
	std::vector<double> compute_wait_suffix() const;
	std::pair<int, int> repair(RemovedVertexStats& removed_stats);
	RoutePerfStats collect_route_stats() const;
};


class Sol :public Tour
{
public:
	Ins* ins;//pointer to instance object
	std::vector<Tour>tours;//solution consist of collection of tours
	std::vector<int>tourindex;//random tour index
	boost::dynamic_bitset<> available;// bitset that states for every vertex if it is still available for inclusion
	int score;//sum of all tour scores
	//methods
	Sol(){}
	~Sol(){}
	Sol(Ins& ins);
	void insert_vertex(Tour &tour, Ins::Vertex* candidate, int position);
	void replace_vertex(Tour &tour, Ins::Vertex* candidate, int position);
	void replace_vertex(Tour &tour, Ins::Vertex* candidate, int position, int breakindex);
	void remove_vertex(Tour &tour, int position);//update score and availability
	void remove_vertices(Tour &tour, int position1, int position2);//update score and availability
	void reset();
	bool check();
	void read_from_file();
	void write_to_file();
	void write_to_cplex(const std::string& dataset);
	void check_availability();
	SolutionPerfStats repair_and_collect_stats();
	friend std::ostream& operator<<(std::ostream& output, Sol& sol);
};