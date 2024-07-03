#pragma once
#include "instance.h"
#include "solution.h"

class Moves
{
public:
	bool insert_nb(Sol& sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	bool exchange_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	void pull_break(Sol& sol, int tour);
	void reschedule_breaks(Sol& sol);
	bool exchange2_nb(Sol& sol);
	bool swap_nb(Sol& sol, int mode);//mode: 0 first improvement, 1 best improvement
	bool two_opt_nb(Sol& sol, int mode);//mode: 0 first improvement, 1 best improvement
	bool move_nb(Sol& sol, int mode=1);//mode: 0 first improvement, 1 best improvement
	bool swap2_nb(Sol & sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	bool one_one_replace_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	bool two_one_replace_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	Ins* ins;
	Moves(Ins& ins) :ins(&ins) {}
};