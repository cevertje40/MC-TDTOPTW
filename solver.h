#pragma once
#include "moves.h"

class Res
{
public:
	Sol sol;
	double time;
	double gap;
	int removed;
	Res() {}
	Res(Sol sol, double& time) :sol(sol), time(time), gap(0.0) {}
	Res(Sol sol, double& time, int& bestknown) :sol(sol), time(time)
	{
		gap = (double(bestknown - sol.score) / bestknown)*100;
	}
	friend std::ostream& operator<<(std::ostream& output, Sol& res);
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
	std::vector<std::vector<double>> tau;
	std::vector<std::vector<double>> eta;
	std::vector<Sol> s;//solution container
	Sol gb;//global best solution
	int iter_nr;//nr of best ant of the iteration
	double iter_score;//score of best ant of the iteration
	void construct(Sol& sol);
public:
	Aco(Ins& ins, double alpha, double beta, double rho, int max_ants, int max_sol, double max_ni_p,double p_best);
	void update_global_best();
	void pheromone_update();
	Res solve(int bestknown=1);
	std::string name = "ACO";
};

class Ils : public Moves
{
private:
	
	int max_it;
	int threshold1;//controls switch to current best when noimpr
	int threshold2;//controls switch to next random route when noimpr
	int threshold3;//controls switch to next random route when noimpr
	Sol s;//current iteration solution
	Sol gb;//global best solution
	void serial_construct(Sol& sol);
	void shake(Sol& sol, int post, int cons);
public:
	Ils(Ins& ins, int max_sol, int threshold1, int threshold2, int threshold3);
	Res solve(int bestknown=1);
	std::string name = "ILS";
};

using Move = std::pair<std::vector<Ins::Vertex*>, std::vector<Ins::Vertex*>>;

class Tabu : public Moves
{
private:
	int max_noimpr;
	int nb_tabu_it;
	void parallel_construct(Sol& sol);
	Sol s;//current iteration solution
	Sol gb;//global best solution
	

public:
	typedef double (Moves::* ScoreFunctionPointer)(double, double, double, double, double, double, double);
	std::vector<ScoreFunctionPointer>ratiofunctions;
	
	Tabu(Ins& ins, int max_noimpr, int nb_tabu_it);
	void perturbe(Sol& sol);
	Res solve(int bestknown = 1);
	std::string name = "Tabu";
};

class HALNS : public Moves
{
	int max_it;//maximum number of iterations
	int T_init;//controls switch to current best when noimpr
	Sol s;//current iteration solution
	Sol gb;//global best solution
	void remove(Sol& sol, int criteria,int amount);//todo
	void insert(Sol& sol, int criteria);//todo
public:
	HALNS(Ins& ins, int max_it, int T_init);
	Res solve(int bestknown = 1);
	std::string name = "HALNS";
};



