#pragma once
#include "instance.h"
#include "solution.h"

class Moves
{
public:
	void best_insert_nb(Sol& sol);
	void pull_break(Sol& sol, int tour);
	void reschedule_breaks(Sol& sol);
	void best_replace_nb(Sol& sol);
	void replace_nb(Sol& sol);
	void exchange(Sol& sol);
	void swap_nb(Sol& sol);
	void two_opt_nb(Sol& sol);
	void best_move_nb(Sol& sol);
	void best_swap2_nb(Sol & sol);
	mt19937 mt;
	Ins* ins;
	Moves(Ins& ins) :ins(&ins) {}
};