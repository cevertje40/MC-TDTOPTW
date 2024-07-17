#pragma once
#include "instance.h"
#include "solution.h"

class One_one_rep_nb
{
public:
	int tour;
	int rempos;
	int inspos;
	Ins::Vertex* inscand;
	double score;
	One_one_rep_nb(int tour, int rempos, int inspos, Ins::Vertex* inscand, double score) :tour(tour), rempos(rempos), inspos(inspos), inscand(inscand), score(score) {}
	friend bool operator< (const One_one_rep_nb& x, const One_one_rep_nb& y)
	{
		return x.score < y.score;
	}
	void execute(Sol& sol) 
	{
		if (rempos != -1)//replace
		{
			sol.removevertex(sol.tours[tour], rempos);
		}
		sol.insertvertex(sol.tours[tour], inscand, inspos);
	}
};

class One_two_rep_nb
{
public:
	int tour;
	int rempos;
	int inspos1;
	int	inspos2;
	Ins::Vertex* inscand1;
	Ins::Vertex* inscand2;
	double score;
	One_two_rep_nb(int tour, int rempos, int inspos1, int inspos2, Ins::Vertex* inscand1, Ins::Vertex* inscand2, double score) :tour(tour), rempos(rempos), inspos1(inspos1), inspos2(inspos2), inscand1(inscand1), inscand2(inscand2), score(score) {}
	friend bool operator< (const One_two_rep_nb& x, const One_two_rep_nb& y)
	{
		return x.score < y.score;
	}
	void execute(Sol& sol)
	{
		sol.removevertex(sol.tours[tour], rempos);
		sol.insertvertex(sol.tours[tour], inscand1, inspos1);
		sol.insertvertex(sol.tours[tour], inscand2, inspos2 + 1);
	}
};

class Two_one_rep_nb
{
public:
	int tour;
	int rempos1;
	int rempos2;
	int inspos;
	Ins::Vertex* inscand;
	double score;
	Two_one_rep_nb(int tour, int rempos1, int rempos2, int inspos, Ins::Vertex* inscand, double score) :tour(tour), rempos1(rempos1), rempos2(rempos2), inspos(inspos), inscand(inscand), score(score) {}
	friend bool operator< (const Two_one_rep_nb& x, const Two_one_rep_nb& y)
	{
		return x.score < y.score;
	}
	void execute(Sol& sol)
	{
		sol.removevertex(sol.tours[tour], rempos1);
		if (rempos1 < rempos2)
		{
			sol.removevertex(sol.tours[tour], rempos2 - 1);
		}
		else
		{
			sol.removevertex(sol.tours[tour], rempos2);
		}
		sol.insertvertex(sol.tours[tour], inscand, inspos);
	}
};

class Moves
{
public:
	bool insert_nb(Sol& sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	bool exchange_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	void pull_break(Sol& sol, int tour);
	void reschedule_breaks(Sol& sol);
	bool exchange2_nb(Sol& sol);
	bool swap_nb(Sol& sol, int mode=1);//mode: 0 first improvement, 1 best improvement
	bool two_opt_nb(Sol& sol, int mode=1);//mode: 0 first improvement, 1 best improvement
	bool move_nb(Sol& sol, int mode=1);//mode: 0 first improvement, 1 best improvement
	bool swap2_nb(Sol& sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	bool one_one_replace_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool two_one_replace_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool one_two_replace(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	Ins* ins;
	//nb generating moves
	boost::heap::priority_queue<One_one_rep_nb> one_one_replace_gen_nb(Sol& sol, int limit);
	boost::heap::priority_queue<One_two_rep_nb> one_two_replace_gen_nb(Sol& sol, int limit);
	boost::heap::priority_queue<Two_one_rep_nb> two_one_replace_gen_nb(Sol& sol, int limit);
	Moves(Ins& ins) :ins(&ins) {}
};

