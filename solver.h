#pragma once
#include "solution.h"

class Solver
{
private:
	void best_insert_nb(Sol& sol);
	void pull_break(Sol& sol, int path);
	void reschedule_breaks(Sol& sol);
	void best_replace_nb(Sol& sol);
	void replace_nb(Sol& sol);
	void exchange(Sol& sol);
	void swap_nb(Sol& sol);
	void two_opt_nb(Sol& sol);
	void ruin_recreate(Sol& sol);
	void relocate_nb(Sol& sol);
public:
	void aco(Ins& instance, double alpha, double beta, double rho,int maxants,int max_it);
};

