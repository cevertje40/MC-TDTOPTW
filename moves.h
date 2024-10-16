#pragma once
#include "instance.h"
#include "solution.h"
#include "tabuvector.h"

class One_one_rep_nb
{
public:
	pair<vector<Ins::Vertex*>,vector<Ins::Vertex*>>move;//out, in
	int tour;
	int rempos;//removal index
	int inspos;//insert index
	Ins::Vertex* inscand;
	double ratio;//ratio increase/decrease
	double score;//score increase/decrease
	One_one_rep_nb(int tour, int rempos, int inspos, Ins::Vertex* inscand, double score, double ratio, pair<vector<Ins::Vertex*>, vector<Ins::Vertex*>>move) :tour(tour), rempos(rempos), inspos(inspos), inscand(inscand), score(score), ratio(ratio), move(move) {}
	friend bool operator< (const One_one_rep_nb& x, const One_one_rep_nb& y)
	{
		return x.ratio < y.ratio;
	}
	void execute(Sol& sol) 
	{
		if (rempos != -1)//replace
		{
			sol.remove_vertex(sol.tours[tour], rempos);
		}
		//otherwise insert move
		sol.insert_vertex(sol.tours[tour], inscand, inspos);
	}
};

class One_two_rep_nb
{
public:
	pair<vector<Ins::Vertex*>, vector<Ins::Vertex*>>move;//out, in
	int tour;
	int rempos;
	int inspos1;
	int	inspos2;
	Ins::Vertex* inscand1;
	Ins::Vertex* inscand2;
	double ratio;//ratio increase/decrease
	double score;//score increase/decrease
	One_two_rep_nb(int tour, int rempos, int inspos1, int inspos2, Ins::Vertex* inscand1, Ins::Vertex* inscand2, double score,double ratio, pair<vector<Ins::Vertex*>, vector<Ins::Vertex*>>move) :tour(tour), rempos(rempos), inspos1(inspos1), inspos2(inspos2), inscand1(inscand1), inscand2(inscand2), score(score), ratio(ratio), move(move) {}
	friend bool operator< (const One_two_rep_nb& x, const One_two_rep_nb& y)
	{
		return x.ratio < y.ratio;
	}
	void execute(Sol& sol)
	{
		sol.remove_vertex(sol.tours[tour], rempos);
		sol.insert_vertex(sol.tours[tour], inscand1, inspos1);
		sol.insert_vertex(sol.tours[tour], inscand2, inspos2);
	}
};

class Two_one_rep_nb
{
public:
	pair<vector<Ins::Vertex*>, vector<Ins::Vertex*>>move;//out, in
	int tour;
	int rempos1;
	int rempos2;
	int inspos;
	Ins::Vertex* inscand;
	double ratio;//ratio increase/decrease
	double score;//score increase/decrease
	Two_one_rep_nb(int tour, int rempos1, int rempos2, int inspos, Ins::Vertex* inscand, double score,double ratio, pair<vector<Ins::Vertex*>, vector<Ins::Vertex*>>move) :tour(tour), rempos1(rempos1), rempos2(rempos2), inspos(inspos), inscand(inscand), score(score), ratio(ratio), move(move) {}
	friend bool operator< (const Two_one_rep_nb& x, const Two_one_rep_nb& y)
	{
		return x.ratio < y.ratio;
	}
	void execute(Sol& sol)
	{
		sol.remove_vertices(sol.tours[tour],rempos1,rempos2);
		sol.insert_vertex(sol.tours[tour], inscand, inspos);
	}
};

class Moves
{
public:
	//local search moves
	bool insert_nb(Sol& sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	bool exchange_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	void pull_break(Sol& sol, int tour);//tries to pull the break forwards
	void reschedule_breaks(Sol& sol);
	bool exchange2_nb(Sol& sol);
	bool swap_nb(Sol& sol, int mode=1);//mode: 0 first improvement, 1 best improvement
	bool two_opt_nb(Sol& sol, int mode=1);//mode: 0 first improvement, 1 best improvement
	bool shift_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool move_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool swap2_nb(Sol& sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	bool one_one_replace_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool two_one_replace_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool one_two_replace(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool or_opt(Sol& sol, int mode = 1);
	Ins* ins;
	//evaluation criteria for nb generators
	double score(double tt, double score, double weight, double volume, double alpha, double beta, double gamma);
	double score_tt(double tt, double score, double weight, double volume, double alpha, double beta, double gamma);
	double score_w(double tt, double score, double weight, double volume, double alpha, double beta, double gammaa);
	double score_v(double tt, double score, double weight, double volume, double alpha, double beta, double gammaa);
	double weighted_ratio(double tt, double score, double weight, double volume, double alpha, double beta, double gamma);
	//nb generators
	boost::heap::priority_queue<One_one_rep_nb> one_one_replace_gen_nb(Sol& sol, TabuVector& tabulist, int globalbest, double alpha, double beta, double gamma, double (Moves::*get_score)(double, double, double, double, double, double, double));
	boost::heap::priority_queue<One_two_rep_nb> one_two_replace_gen_nb(Sol& sol, TabuVector& tabulist, int globalbest, double alpha, double beta, double gamma, double (Moves::*get_score)(double, double, double, double, double, double, double));
	boost::heap::priority_queue<Two_one_rep_nb> two_one_replace_gen_nb(Sol& sol, TabuVector& tabulist, int globalbest, double alpha, double beta, double gamma, double (Moves::*get_score)(double, double, double, double, double, double, double));
	boost::heap::priority_queue<One_one_rep_nb> one_one_replace_gen_nb(Sol& sol, TabuVector& tabulist,int globalbest);
	boost::heap::priority_queue<One_two_rep_nb> one_two_replace_gen_nb(Sol& sol, TabuVector& tabulist,int globalbest);
	boost::heap::priority_queue<Two_one_rep_nb> two_one_replace_gen_nb(Sol& sol, TabuVector& tabulist,int globalbest);
	//constructor
	Moves(Ins& ins) :ins(&ins) {}
};

