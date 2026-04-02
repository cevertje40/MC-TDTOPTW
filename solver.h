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

enum class SelectionOp
{
	Random = 1,
	HighestScore = 2,
	HighestScoreToBurden = 3,
	LowestResource = 4,
	RegretTime = 5
};

enum class InsertionOp
{
	BestPosition = 1,
	FirstFeasible = 2,
	LastFeasible = 3,
	RandomFeasible = 4,
	LeastLoadedBestPosition = 5
};

template <size_t N>
int weighted_pick(const std::array<double, N>& w, std::mt19937& engine)
{
	double sum = 0.0;
	for (double x : w) sum += x;
	if (sum <= 0.0)
	{
		std::uniform_int_distribution<int> dist(0, (int)N - 1);
		return dist(engine);
	}

	double u = std::uniform_real_distribution<double>(0.0, sum)(engine);
	double acc = 0.0;
	for (int i = 0; i < (int)N; ++i)
	{
		acc += w[i];
		if (u <= acc) return i;
	}
	return (int)N - 1;
}

class Alns : public Moves
{
private:
	int max_it;
	Sol s;//current iteration solution
	Sol gb;//global best solution

	std::array<double, 7> rem_w;
	std::array<double, 5> sel_w;
	std::array<double, 5> ins_w;

	std::array<double, 7> rem_score;
	std::array<double, 5> sel_score;
	std::array<double, 5> ins_score;

	std::array<int, 7> rem_used;
	std::array<int, 5> sel_used;
	std::array<int, 5> ins_used;
	
	std::vector<std::pair<int, int>> collect_removable_positions(Sol& sol);
	const RemovedCustomer* find_removed_info(Ins::Vertex* y, const std::vector<RemovedCustomer>& removed) const;
	//removal operators
	std::vector<RemovedCustomer> random_remove(Sol& sol, int beta);
	std::vector<RemovedCustomer> worst_remove_burden(Sol& sol, int beta);
	std::vector<RemovedCustomer> largest_time_saving_remove(Sol& sol, int beta);
	std::vector<RemovedCustomer> sequence_remove(Sol& sol, int beta);
	std::vector<RemovedCustomer> break_neighborhood_remove(Sol& sol, int beta);
	std::vector<RemovedCustomer> largest_demand_remove(Sol& sol, int beta);
	std::vector<RemovedCustomer> lowest_profit_remove(Sol& sol, int beta);
	std::vector<RemovedCustomer> largest_service_time_remove(Sol& sol, int beta);
	std::vector<RemovedCustomer> random_route_remove(Sol& sol, int beta);


	
	//selection operators
	std::vector<Ins::Vertex*> collect_available_customers(const Sol& sol) const;

	Ins::Vertex* random_selection_from_pool(const std::vector<Ins::Vertex*>& pool) const;
	Ins::Vertex* highest_score_selection_from_pool(const std::vector<Ins::Vertex*>& pool) const;
	Ins::Vertex* highest_score_to_burden_selection_from_pool(const std::vector<Ins::Vertex*>& pool) const;
	Ins::Vertex* lowest_resource_selection_from_pool(const std::vector<Ins::Vertex*>& pool) const;
	Ins::Vertex* dynamic_travel_time_profit_selection_from_pool(const Sol& sol, const std::vector<Ins::Vertex*>& pool, const std::vector<RemovedCustomer>& removed) const;

	Ins::Vertex* apply_selection_operator_from_pool(const Sol& sol, const std::vector<Ins::Vertex*>& pool, SelectionOp sel_op, const std::vector<RemovedCustomer>& removed) const;

	//insertion operators
	std::vector<double> collect_feasible_insertion_shifts(const Sol& sol, Ins::Vertex* y, const std::vector<RemovedCustomer>& removed) const;
	bool evaluate_insertion_position(const Sol& sol, Ins::Vertex* y, int d, int j, const std::vector<RemovedCustomer>& removed, GreedyInsertion& out) const;

	GreedyInsertion best_position_insertion(const Sol& sol, Ins::Vertex* y, const std::vector<RemovedCustomer>& removed) const;
	GreedyInsertion first_feasible_insertion(const Sol& sol, Ins::Vertex* y, const std::vector<RemovedCustomer>& removed) const;
	GreedyInsertion last_feasible_insertion(const Sol& sol, Ins::Vertex* y, const std::vector<RemovedCustomer>& removed) const;
	GreedyInsertion random_feasible_insertion(const Sol& sol, Ins::Vertex* y, const std::vector<RemovedCustomer>& removed) const;
	GreedyInsertion least_loaded_best_position_insertion(const Sol& sol, Ins::Vertex* y, const std::vector<RemovedCustomer>& removed) const;

	GreedyInsertion apply_insertion_operator(const Sol& sol, Ins::Vertex* y, InsertionOp op, const std::vector<RemovedCustomer>& removed) const;//dispatcher for insertion operators

	void apply_local_search(Sol& sol);

	void repair_with_selection_and_insertion(Sol& sol, SelectionOp sel_op, InsertionOp ins_op, const std::vector<RemovedCustomer>& removed);



	//acceptance functions
	bool accept_candidate(const Sol& cur, const Sol& cand, double T);

	void update_operator_weights();


public:
	Alns(Ins& ins, int max_it);
	Res solve(int bestknown = 1);

};


