#pragma once
#include "moves.h"
#include "elitepool.h"


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
	double EMA_RHO; //	learning rate for exponential moving average of gains
	double THRESH;  // start biasing after 80% utilization
	double GAMMA;  // bias strength
	
	Sol s;//current iteration solution
	Sol gb;//global best solution
	

public:
	
	Tabu(Ins& ins, int max_noimpr, int nb_tabu_it,double ema_rho, double thresh, double gamma);
	void perturbe(Sol& sol);
	Res solve(int bestknown = 1);
	std::string name = "Tabu";
};

struct RemovedCustomer
{
	Ins::Vertex* v = nullptr;
	int oldTour = -1;
	int oldPos = -1;
	int oldPred = -1;
	int oldSucc = -1;
};

struct GreedyInsertion
{
	bool feasible = false;
	int vertex_idx = -1;
	int tour_idx = -1;
	int position = -1;   // insert after seq[position]
	double key = -DBL_MAX;
	double shift = 0.0;
};

class Alns : public Moves
{
private:
	int max_it;
	Sol s;//current iteration solution
	Sol gb;//global best solution
	
	std::vector<std::pair<int, int>> collect_removable_positions(Sol& sol);
	std::vector<RemovedCustomer> random_remove_1(Sol& sol);
	std::vector<RemovedCustomer> random_remove_2(Sol& sol);
	std::vector<RemovedCustomer> worst_remove_1(Sol& sol);
	std::vector<RemovedCustomer> worst_remove_1_resource_time(Sol& sol);
	GreedyInsertion best_insertion_for_vertex(const Sol& sol, Ins::Vertex* y, const std::vector<RemovedCustomer>& removed);
	std::vector<Ins::Vertex*> collect_available_customers(const Sol& sol);
	const RemovedCustomer* find_removed_info(Ins::Vertex* y,const std::vector<RemovedCustomer>& removed) const;
	void greedy_repair(Sol& sol, const std::vector<RemovedCustomer>& removed);
	bool accept_candidate(const Sol& cur, const Sol& cand, double T);
public:
	Alns(Ins& ins, int max_it);
	Res solve(int bestknown = 1);

};


