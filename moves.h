#pragma once
#include "instance.h"
#include "solution.h"

class Moves
{
public:
	void insert_nb(Sol& sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	void replace_nb(Sol& sol, int mode = 1);//mode: 0 first improvement, 1 best improvement
	void pull_break(Sol& sol, int tour);
	void reschedule_breaks(Sol& sol);
	void exchange(Sol& sol);
	void swap_nb(Sol& sol, int mode);//mode: 0 first improvement, 1 best improvement
	void two_opt_nb(Sol& sol, int mode);//mode: 0 first improvement, 1 best improvement
	void move_nb(Sol& sol, int mode=1);//mode: 0 first improvement, 1 best improvement
	void swap2_nb(Sol & sol,int mode=1);//mode: 0 first improvement, 1 best improvement
	mt19937 mt;
	Ins* ins;
	Moves(Ins& ins) :ins(&ins) {}
};