#pragma once
#include "moves.h"

class Res
{
public:
	int score;
	double time;
	int bestknown;
	double gap;
	Res() {}
	Res(int& score, double& time) :score(score), time(time), bestknown(0), gap(0.0) {}
	Res(int& score, double& time,int& bestknown):score(score), time(time), bestknown(bestknown)
	{
		gap = (double(bestknown - score) / bestknown)*100;
	}
	friend ostream& operator<<(ostream& output, Sol& res);
};

class Aco: public Moves
{
private:
	double alpha;
	double beta;
	double rho;
	int max_ants;
	int max_it;
	int ni;//number of non improvement iterations
	int max_ni;//number of iterations allowed without improvement before pheromone reset
	double p_best;//controls pheromone update process
	vector<vector<double>> tau;
	vector<vector<double>> eta;
	vector<Sol> s;//solution container
	Sol gb;//global best solution
	int iter_nr;//nr of best ant of the iteration
	double iter_score;//score of best ant of the iteration
	void construct(Sol& sol);
public:
	Aco(Ins& ins, double alpha, double beta, double rho, int max_ants, int max_sol, double max_ni_p,double p_best);
	void update_global_best();
	void pheromone_update();
	Res solve(int bestknown=1);
};

class Ils : public Moves
{
private:
	
	int max_it;
	int ni;//number of non improvement iterations
	int threshold1;//controls switch to current best when noimpr
	int threshold2;//controls switch to next random route when noimpr
	int threshold3;//controls switch to next random route when noimpr
	Sol s;//current iteration solution
	Sol gb;//global best solution
	int iter_nr;//nr of best ant of the iteration
	double iter_score;//score of best ant of the iteration
	void serial_construct(Sol& sol);
	void parallel_construct(Sol& sol);
	void shake(Sol& sol, int post, int cons);
public:
	Ils(Ins& ins, int max_sol, int threshold1, int threshold2, int threshold3);
	Res solve(int bestknown=1);
};


