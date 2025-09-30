#pragma once
#include "instance.h"

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
	std::pair<int,int> repair();//repairs solution by removed last regular vertex, return score decrease
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
	void write_to_cplex();
	void check_availability();
	int repair();
	friend std::ostream& operator<<(std::ostream& output, Sol& sol);
};