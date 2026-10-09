#include "solver.h"

using namespace std;

Aco::Aco(Ins& ins, double alpha, double beta, double rho, int max_ants, int max_sol, double max_ni_p, double p_best): Moves(ins), alpha(alpha), beta(beta), rho(rho), max_ants(max_ants), ni(0), p_best(p_best)
{
	max_it = int(max_sol / max_ants);
	max_ni = int(max_ni_p * max_it);
	//creation of eta and theta
	tau.resize(ins.maxvertices - 1);
	eta.resize(ins.maxvertices - 1);
	for (int i = 0; i < ins.maxvertices - 1; ++i)
	{
		tau[i].resize(ins.maxvertices);
		eta[i].resize(ins.maxvertices);
		for (int j = 0; j < ins.maxvertices; ++j)
		{
			tau[i][j] = 1.0;
			if (i != j)
			{
				//double consumption = ins.v[i].con[j]->determin / ins.t[0].T_max;
				double consumption = (((ins.v[i].con[j]->determin+ins.v[j].serv) / ins.t[0].T_max) + (ins.v[j].weight / ins.t[0].W_max) + (ins.v[j].volume / ins.t[0].V_max)) / 3;
				eta[i][j] = std::max(1,ins.v[j].score) / consumption;
			}
			else
				eta[i][j] = 0.0;
		}
	}
	//creation of solutions
	gb = Sol(ins);//ant class object
	s.resize(max_ants);//array of ant class objects
	#pragma omp parallel for// parallel first touch to increase speed
	for (int i = 0; i < max_ants; ++i) 
	{
		s[i] = Sol(ins);
	}// end for all ants
	iter_nr = 0;
	iter_score = -1;
}

void Aco::update_global_best()
{
	for (int ant = 0; ant < max_ants; ++ant) //find iteration best
	{
		if (iter_score < s[ant].score)
		{
			iter_score = s[ant].score;
			iter_nr = ant;
		}
	}
	if (gb.score < iter_score)//update best solution
	{
		gb = s[iter_nr];
		ni = 0;
	}
	else//increase number of iterations without improvement
	{
		++ni;
	}
}

void Aco::pheromone_update()
{
	double tau_max = (gb.score / ins->maxscore) / (1 - rho);
	double tau_min = (1 - (pow(p_best, 1.0 / ins->maxvertices))) / (((double(ins->maxvertices) / 2) - 1) * (pow(p_best, 1.0 / ins->maxvertices))) * tau_max;
	if (ni < max_ni)
	{
		//augment edges to iteration's best solution
		double quality = s[iter_nr].score / ins->maxscore;
		//decrease and check limits
		for (int t = 0; t < ins->maxtours; ++t)
		{
			int end = (int)s[iter_nr].tours[t].seq.size() - 1;
			for (int i = 0; i < end; ++i)
			{
				tau[s[iter_nr].tours[t].seq[i]->index][s[iter_nr].tours[t].seq[i + 1]->index] += quality;
			}
		}
		for (int i = 0; i < ins->maxvertices - 1; ++i)
		{
			for (int j = 0; j < ins->maxvertices; ++j)
			{
				tau[i][j] *= (1 - rho);
				if (tau[i][j] < tau_min)
					tau[i][j] = tau_min;
				else if (tau[i][j] > tau_max)
					tau[i][j] = tau_max;
			}
		}
	}
	else//if number of non improving iterations is achieved reset pheromones
	{
		ni = 0;
		//reset all edges
		for (int i = 0; i < ins->maxvertices - 1; ++i)
		{
			for (int j = 0; j < ins->maxvertices; ++j)
			{
				tau[i][j] = tau_max;
			}
		}
	}
	iter_nr = 0;
	iter_score = 0;
}

void Aco::construct(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int t = sol.tourindex[d];
		Sol::Tour& tour = sol.tours[sol.tourindex[d]];
		bool breaktaken = false;
		while (tour.seq.back()->index != ins->maxvertices - 1) //until one solution is full=>sequential procedure
		{
			vector<double> prob_v;
			vector<double> temp_traveltime;
			vector<int> temp_action;
			Ins::Vertex* last = tour.seq.back();
			prob_v.resize(ins->maxvertices);
			temp_traveltime.resize(ins->maxvertices);
			temp_action.resize(ins->maxvertices);
			for (int i = 0; i < ins->maxvertices; ++i)
			{//check 3 constraints for potential vertex to include: availability,weight,volume
				if (i != ins->maxvertices - 1)
				{
					if (sol.available[i])
					{
						if ((tour.volume + ins->v[i].volume <= ins->t[t].V_max) && (tour.weight + ins->v[i].weight <= ins->t[t].W_max))
						{
							prob_v[i] = 1;
						}
					}
					else
					{
						prob_v[i] = 0;
					}
				}
				else
				{// always add enddepot als potential vertex to include
					prob_v[i] = 1;
				}
			}
			//travel time check for potential inclusions
			bool enddepot = true;
			for (int i = 0; i < ins->maxvertices; ++i)
			{
				if (prob_v[i] != 0)
				{
					double currenttime = tour.deptime.back() + ins->t[t].EDT;
					Ins::Vertex* neighbor = &ins->v[i];
					double arrivaltime = ins->arrival_time(last->con[neighbor->index], currenttime);
					int brk = -1;
					//als je na breakstart aankomt of moet wachten 
					if ((breaktaken == false) && ((arrivaltime >= ins->breakstart)||(neighbor->LTW[t]-ins->breakdur>=ins->breakstart)))
					{
						arrivaltime += ins->breakdur;
						brk = 1;
					}
					else
					{
						brk = 0;
					}
					if (arrivaltime < neighbor->LTW[t])
					{
						prob_v[i] = 1 - (double(neighbor->LTW[t] - arrivaltime) / ins->t[t].T_max);//wachttijd meegeven
						arrivaltime = neighbor->LTW[t];//wachten als je te vroeg bent	
					}
					if (arrivaltime > neighbor->UTW[t])
					{
						prob_v[i] = 0;
						continue;//mag je al stoppen met rekenen
					}
					arrivaltime += neighbor->serv;
					double enddepottime = ins->arrival_time(neighbor->con[ins->maxvertices - 1], arrivaltime);
					if ((breaktaken == false) && (brk == 0))//no break taken yet and you are not going to break this iter
					{
						enddepottime += ins->breakdur;//take break at end depot
					}
					if (enddepottime > ins->t[t].LAT)//enddepot heeft geen service time
					{
						prob_v[i] = 0;//discard vertex if infeasible
					}
					else
					{
						if (neighbor->index != ins->maxvertices - 1)
						{
							enddepot = false;//if another point than the enddepot is found
						}
						temp_traveltime[i] = arrivaltime - ins->t[t].EDT;
						temp_action[i] = brk;
					}
				}//end if prob
			}//end for
			//don't include enddepot if there are still feasible vertices
			if (enddepot == false)
			{
				prob_v[ins->maxvertices - 1] = 0;
			}
			// probability calculation
			double denominator = 0.0;
			for (int i = 0; i < ins->maxvertices; ++i)
			{
				if (prob_v[i] != 0)
				{
					double temp1 = tau[last->index][i];
					double temp2 = eta[last->index][i];
					for (int a = 0; a < alpha; ++a)
					{
						prob_v[i] *= temp1;
					}
					for (int b = 0; b < beta; ++b)
					{
						prob_v[i] *= temp2;
					}
					denominator += prob_v[i];
				}
			}
			for (int i = 0; i < ins->maxvertices; ++i)
			{
				if (prob_v[i] != 0)
					prob_v[i] /= denominator;
			}
			//selection procedure based on prob_v and random number
			if (enddepot == true)
			{
				prob_v[ins->maxvertices - 1] = 2;
			}
			double total = 0.0;
			int sel = 0;
			uniform_real_distribution<double> dist(0, 1);
			double random = dist(engine);
			while (random >= total)
			{
				total += prob_v[sel];
				++sel;
			}
			//cout << "hier" << endl;
			//cin >> sel;
			//sel += 1;
			tour.seq.push_back(&ins->v[sel - 1]);
			tour.deptime.push_back(temp_traveltime[sel - 1]);
			sol.available[sel - 1] = false;
			tour.score += ins->v[sel - 1].score;
			sol.score += ins->v[sel - 1].score;
			tour.max_shift.push_back(0);//dummy die dan in calc max shift upgedate wordt
			tour.volume += ins->v[sel - 1].volume;
			tour.weight += ins->v[sel - 1].weight;
			tour.action.push_back(temp_action[sel - 1]);
			if (temp_action[sel - 1] == 1)
			{
				breaktaken = true;
				tour.breakindex = int(tour.seq.size()) - 1;
			}
			prob_v[sel - 1] = 0;
			//++counter;
		}// end while

		sol.available[ins->maxvertices - 1] = true;// end depot moet terug available zijn voor de volgende tour
		if (breaktaken == false)
		{
			tour.action.back() = 1;
			tour.deptime.back() += ins->breakdur;
			tour.breakindex = int(tour.seq.size()) - 1;
		}
		tour.calc_maxshift();
	}//end for all d
	//make end depot unaivailable for other moves
	sol.available[ins->maxvertices - 1] = false;
	//sol.check();
}

ostream& operator<<(ostream& output, Res& res)
{
	output << "score: " << res.sol.score << " after: " << res.time << " gap: " << res.gap << "\n";
	return output;
}

Res Aco::solve(int bestknown)
{
	clock_t start, end;
	start = clock();
	for (int iter = 0; iter < max_it; ++iter)//iteration loop
	{
		#pragma omp parallel
		{//start parallel session
			#pragma omp for nowait
			for (int ant = 0; ant < max_ants; ++ant) // for all ants
			{
				s[ant].reset();
				construct(s[ant]);
				//two_opt_nb(s[ant], 1);
				swap_nb(s[ant],1);
				swap2_nb(s[ant], 1);
				relocate_nb(s[ant], 1);
				insert_nb(s[ant], 1);
				one_one_replace_nb(s[ant], 1);
				//two_one_replace_nb(s[ant], 1);
				//one_two_replace(s[ant], 1);
				//exchange_nb(s[ant], 1);
				//two_opt_nb(s[ant], 1);
				//swap_nb(s[ant], 1);
				//swap2_nb(s[ant], 1);
				//move_nb(s[ant], 1);
				//insert_nb(s[ant], 1);
				//exchange_nb(s[ant], 1);
				//insert_nb(s[ant], 1);
			}
		}
		//cout << iter << endl;
		update_global_best();//best solution is stored
		pheromone_update();// iteration best solution its arcs are augmented
	}
	end = clock();
	double cpuTime;
	cpuTime = difftime(end, start) / CLOCKS_PER_SEC;
	gb.check();
	std::cout << gb << endl;
	//cout << "best score: "<<gb.score<<"after: "<<cpuTime << endl;
	return Res(gb, cpuTime, bestknown);
}

void Ils::serial_construct(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int t = sol.tourindex[d];
		Sol::Tour& tour = sol.tours[sol.tourindex[d]];
		bool breaktaken = false;
		while (tour.seq.back()->index != ins->maxvertices - 1) //until one solution is full=>sequential procedure
		{
			vector<double> prob_v;
			vector<double> temp_traveltime;
			vector<int> temp_action;
			Ins::Vertex* last = tour.seq.back();
			prob_v.resize(ins->maxvertices);
			temp_traveltime.resize(ins->maxvertices);
			temp_action.resize(ins->maxvertices);
			for (int i = 0; i < ins->maxvertices; ++i)
			{//check 3 constraints for potential vertex to include: availability,weight,volume
				if (i != ins->maxvertices - 1)
				{
					if (sol.available[i])
					{
						if ((tour.volume + ins->v[i].volume <= ins->t[t].V_max) && (tour.weight + ins->v[i].weight <= ins->t[t].W_max))
						{
							prob_v[i] = 1;
						}
					}
					else
					{
						prob_v[i] = 0;
					}
				}
				else
				{// always add enddepot als potential vertex to include
					prob_v[i] = 1;
				}
			}
			//travel time check for potential inclusions
			bool enddepot = true;
			for (int i = 0; i < ins->maxvertices; ++i)
			{
				if (prob_v[i] != 0)
				{
					double currenttime = tour.deptime.back() + ins->t[t].EDT;
					Ins::Vertex* neighbor = &ins->v[i];
					double arrivaltime = ins->arrival_time(last->con[neighbor->index], currenttime);
					int brk = -1;
					//als je na breakstart aankomt of moet wachten 
					if ((breaktaken == false) && ((arrivaltime >= ins->breakstart) || (neighbor->LTW[t] - ins->breakdur >= ins->breakstart)))
					{
						arrivaltime += ins->breakdur;
						brk = 1;
					}
					else
					{
						brk = 0;
					}
					if (arrivaltime < neighbor->LTW[t])
					{
						prob_v[i] = 1 - (double(neighbor->LTW[t] - arrivaltime) / ins->t[t].T_max);//wachttijd meegeven
						arrivaltime = neighbor->LTW[t];//wachten als je te vroeg bent	
					}
					if (arrivaltime > neighbor->UTW[t])
					{
						prob_v[i] = 0;
						continue;//mag je al stoppen met rekenen
					}
					arrivaltime += neighbor->serv;
					double enddepottime = ins->arrival_time(neighbor->con[ins->maxvertices - 1], arrivaltime);
					if ((breaktaken == false) && (brk == 0))//no break taken yet and you are not going to break this iter
					{
						enddepottime += ins->breakdur;//take break at end depot
					}
					if (enddepottime > ins->t[t].LAT)//enddepot heeft geen service time
					{
						prob_v[i] = 0;//discard vertex if infeasible
					}
					else
					{
						if (neighbor->index != ins->maxvertices - 1)
						{
							enddepot = false;//if another point than the enddepot is found
						}
						temp_traveltime[i] = arrivaltime - ins->t[t].EDT;
						temp_action[i] = brk;
					}
				}//end if prob
			}//end for
			//don't include enddepot if there are still feasible vertices
			if (enddepot == false)
			{
				prob_v[ins->maxvertices - 1] = 0;
			}
			// probability calculation
			double denominator = 0.0;
			for (int i = 0; i < ins->maxvertices; ++i)
			{
				if (prob_v[i] != 0)
				{
					double tt = (last->con[i]->determin) / ins->t[tour.index].T_max;
					double w = (ins->v[i].weight) / ins->t[tour.index].W_max;
					double vol = (ins->v[i].volume) / ins->t[tour.index].V_max;
					double consumption = (tt + w + vol) / 3.0;
					double score_pos = std::max(1.0, static_cast<double>(ins->v[i].score));
					// Avoid divide-by-zero
					constexpr double eps = 1e-12;
					prob_v[i] = score_pos / std::max(consumption, eps);
					denominator += prob_v[i];
				}
			}
			for (int i = 0; i < ins->maxvertices; ++i)
			{
				if (prob_v[i] != 0)
					prob_v[i] /= denominator;
			}
			//selection procedure based on prob_v and random number
			if (enddepot == true)
			{
				prob_v[ins->maxvertices - 1] = 2;
			}
			double total = 0.0;
			int sel = 0;
			uniform_real_distribution<double> dist(0, 1);
			double random = dist(engine);
			while (random >= total)
			{
				total += prob_v[sel];
				++sel;
			}
			//cout << "hier" << endl;
			//cin >> sel;
			//sel += 1;
			tour.seq.push_back(&ins->v[sel - 1]);
			tour.deptime.push_back(temp_traveltime[sel - 1]);
			sol.available[sel - 1] = false;
			tour.score += ins->v[sel - 1].score;
			sol.score += ins->v[sel - 1].score;
			tour.max_shift.push_back(0);//dummy die dan in calc max shift upgedate wordt
			tour.volume += ins->v[sel - 1].volume;
			tour.weight += ins->v[sel - 1].weight;
			tour.action.push_back(temp_action[sel - 1]);
			if (temp_action[sel - 1] == 1)
			{
				breaktaken = true;
				tour.breakindex = int(tour.seq.size()) - 1;
			}
			prob_v[sel - 1] = 0;
			//++counter;
		}// end while

		sol.available[ins->maxvertices - 1] = true;// end depot moet terug available zijn voor de volgende tour
		if (breaktaken == false)
		{
			tour.action.back() = 1;
			tour.deptime.back() += ins->breakdur;
			tour.breakindex = int(tour.seq.size()) - 1;
		}
		tour.calc_maxshift();
	}//end for all d
	//make end depot unaivailable for other moves
	sol.available[ins->maxvertices - 1] = false;
	//sol.check();
}//end serial construct

void Ils::shake(Sol& sol, int post, int cons)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int t = sol.tourindex[d];
		Sol::Tour& tour = sol.tours[sol.tourindex[d]];
		if (tour.seq.size() > 2)//only shake tours containing regular vertices
		{
			cons = min(cons, (int)tour.seq.size() - 2);//const is only capped by maxpaths not the maximum size of the route under consideration
			bool breakupdate = false;
			if (post <= tour.breakindex)
			{
				breakupdate = true;
			}
			for (int i = 0; i < cons; ++i)
			{	
				//removal has reached last regular vertex, continue from beginning
				if (post >= tour.seq.size() - 1)
				{
					post = 1;
					breakupdate = true;//if you remove the first vertex of the route a breakupdate is required
				}
				//remove vertices one by one at position post
				sol.score -= tour.seq[post]->score;
				tour.score -= tour.seq[post]->score;
				tour.weight -= tour.seq[post]->weight;
				tour.volume -= tour.seq[post]->volume;
				sol.available[tour.seq[post]->index] = true;
				tour.seq.erase(tour.seq.begin() + post);
				tour.deptime.erase(tour.deptime.begin() + post);
				tour.max_shift.erase(tour.max_shift.begin() + post);
				tour.action.erase(tour.action.begin() + post);
			}
			if (breakupdate)
			{
				tour.update_break();
			}
			else
			{
				tour.update(post - 1, int(tour.seq.size()));
			}
			//tour.calc_maxshift();
		}//end if tour is not empty
	}//end for all tours
}//end shake

Ils::Ils(Ins& ins, int max_sol, int threshold1, int threshold2, int threshold3): Moves(ins), threshold1(threshold1), threshold2(threshold2), threshold3(threshold3)
{
	max_it = max_sol;
	gb = Sol(ins);//best sol
	s = Sol(ins);//iter sol
}

Res Ils::solve(int bestknown)
{
	clock_t start, end;
	start = clock();
	s.reset();
	//parallel_construct(s);
	serial_construct(s);
	int noimpr = 0;
	int post = 1;//post -->position to start the removal
	int cons = 1;//cons -->amount of vertices to be removed
	for (int iter = 0; iter < max_it; ++iter)//iteration loop
	{
		//cout << iter << endl;
		if ((noimpr > threshold2) && ((noimpr + 1) % threshold3 == 0))
		{
			next_permutation(s.tourindex.begin(), s.tourindex.end());//permutates the path index that is stored in the solution
		}
		else
		{
			shake(s, post, cons);
		}
		two_opt_nb(s,1);
		swap_nb(s,1);
		swap2_nb(s,1);
		relocate_nb(s,1);
		insert_nb(s,1);
		one_one_replace_nb(s, 1);
		//two_one_replace_nb(s, 1);
		//one_two_replace(s, 1);
		if (s.score > gb.score)
		{
			gb = s;
			cons = 1;
			noimpr = 0;
		}
		else
		{
			++noimpr;
		}
		if ((noimpr + 1) % threshold1 == 0)
		{
			s = gb;
			cons = 1;
		}
		post = post + cons;
		if (iter % 2 == 0)
		{
			++cons;
		}
		int maxpathsize = 0;
		int minpathsize = INT_MAX;
		for (int t = 0; t < ins->maxtours; ++t)
		{
			maxpathsize = max(maxpathsize, (int)s.tours[t].seq.size() - 2);
			if (s.tours[t].seq.size() > 2)//only routes that contain regular vertices are checked for minpath size
			{
				minpathsize = min(minpathsize, (int)s.tours[t].seq.size() - 2);
			}
		}
		if (cons > maxpathsize)
		{
			cons = 1;
		}
		while (post > minpathsize)
		{
			post -=minpathsize;
		}
	}
	end = clock();
	double cpuTime;
	cpuTime = difftime(end, start) / CLOCKS_PER_SEC;
	gb.check();
	std::cout << gb << endl;
	//cout << "best score: "<<gb.score<<"after: "<<cpuTime << endl;
	return Res(gb, cpuTime, bestknown);
}


Tabu::Tabu(Ins& ins, int max_noimpr, int nb_tabu_it, double ema_rho, double thresh, double gamma) : Moves(ins), max_noimpr(max_noimpr), nb_tabu_it(nb_tabu_it), ema_rho(ema_rho), thresh(thresh), gamma(gamma)
{
	gb = Sol(ins);//best sol
	s = Sol(ins);//iter sol
}

void Tabu::perturbe(Sol& sol)
{
	constexpr double removalRate = 0.02; // Remove up to 20% of vertices from a tour
	std::uniform_real_distribution<double> prob(0.0, 1.0);
	std::uniform_int_distribution<> offset(1, 2); // how far from the worst to go
	for (int d = 0; d < ins->maxtours; ++d)
	{
		Sol::Tour& tour = sol.tours[d];
		int end = (int)tour.seq.size();
		// Skip if tour is tiny
		if (end <= 3) continue;
		// Identify vertex inefficiencies
		struct Entry {
			int pos;
			double inefficiency;
		};
		std::vector<Entry> ineff;
		for (int i = 1; i < end - 1; ++i)
		{
			Ins::Vertex* a = tour.seq[i - 1];
			Ins::Vertex* b = tour.seq[i];
			double usedTT = (tour.deptime[i] - (b->serv + tour.action[i] * ins->breakdur)) - tour.deptime[i - 1];
			double minTT = a->con[b->index]->determin;
			double ratio = usedTT / max(1.0, minTT);
			ineff.push_back({ i, ratio });
		}
		// Sort by inefficiency (descending)
		std::sort(ineff.begin(), ineff.end(), [](const Entry& lhs, const Entry& rhs) {
			return lhs.inefficiency > rhs.inefficiency;
			});
		int numToRemove = max(1, int(removalRate * ineff.size()));
		unordered_set<int> removed;
		for (int i = 0; i < numToRemove; ++i)
		{
			int index = ineff[i + offset(engine) % min(3, (int)ineff.size())].pos;
			// Avoid removing the same twice
			if (removed.count(index)) continue;
			sol.remove_vertex(tour, index);
			removed.insert(index);
		}
	}
}

template<typename MoveType>
void executeMove(boost::heap::priority_queue<MoveType>& nb, Sol& s, TabuVector& tabulist) 
{
	auto exec_nb = nb.top();
	//cout << nb.top().score<<" " << nb.top().tour << " " << nb.top().ratio << endl;;
	exec_nb.execute(s);
	for (int i = 0; i < exec_nb.move.first.size(); ++i)
	{
		tabulist.addTabu(exec_nb.move.first[i]->index,exec_nb.tour);//out
	}
	for (int i = 0; i < exec_nb.move.second.size(); ++i)
	{
		tabulist.addTabu(exec_nb.move.second[i]->index,exec_nb.tour);//in
	}
}

struct Pressures { double time_p = 0, weight_p = 0, volume_p = 0; };

inline Pressures compute_pressures(const Sol& s) {
	double t_p = 0.0, w_p = 0.0, v_p = 0.0;

	for (size_t d = 0; d < s.tours.size(); ++d) {
		const Tour& tour = s.tours[d];

		// Safety: empty or degenerate tour
		const double used_time = tour.deptime.back();  // hours since EDT
		const double used_w = std::max(0.0, tour.weight);                           // tons
		const double used_v = std::max(0.0, tour.volume);                           // m^3

		const double T = std::max(1e-9, s.ins->t[d].T_max);
		const double W = std::max(1e-9, s.ins->t[d].W_max);
		const double V = std::max(1e-9, s.ins->t[d].V_max);

		t_p = std::max(t_p, std::min(1.0, used_time / T));
		w_p = std::max(w_p, std::min(1.0, used_w / W));
		v_p = std::max(v_p, std::min(1.0, used_v / V));
	}
	return { t_p, w_p, v_p };
}

Res Tabu::solve(int bestknown, double max_time_sec, int forced_crit, int nb_mask)
{
	using clock = std::chrono::steady_clock;
	auto t0 = clock::now();
	s.reset();

	auto elapsed_sec = [&]() -> double
		{
			return std::chrono::duration<double>(clock::now() - t0).count();
		};

	auto time_up = [&]() -> bool
		{
			return elapsed_sec() >= max_time_sec;
		};

	parallel_construct(s);
	//s.read_from_file("20.3.1.3.txt");
	//s.check();
	gb = s;

	ElitePool elites;
	elites.consider(s, ins->maxvertices);

	// criteria selector based on constraint pressure
	enum Crit { SCORE = 0, TIME = 1, VOLUME = 2, WEIGHT = 3, N_CRIT = 4 };
	struct CritStats { double ema_gain = 0.0; long used = 0; };
	std::array<CritStats, N_CRIT> crit_stats{};

	const double EPS_GAIN = 1e-6; // to avoid zero gains

	TabuVector tabulist(ins->maxtours, ins->maxvertices, nb_tabu_it);

	int noimpr = 0;
	int iter = 0;
	int nonb1 = 0;
	int nonb2 = 0;
	int nonb3 = 0;
	int no_feasible_moves_in_a_row = 0;
	int restart_count = 0;
	double restart_dist_sum = 0.0;

	// enabled replacement neighborhoods
	std::vector<int> enabled_nb;
	if (nb_mask & 1) enabled_nb.push_back(1); // 1-1
	if (nb_mask & 2) enabled_nb.push_back(2); // 2-1
	if (nb_mask & 4) enabled_nb.push_back(3); // 1-2

	if (enabled_nb.empty())
	{
		throw std::runtime_error("Tabu::solve called with nb_mask disabling all replacement neighborhoods.");
	}

	while (noimpr < max_noimpr && !time_up())
	{
		if (gb.score == ins->maxscore)
		{
			break;
		}

		Sol backup = s;
		int prev_score = s.score;
		bool moved = false;

		// --- criterion selection ---
		auto P = compute_pressures(s);

		auto ctx_mult = [&](int crit) -> double
			{
				if (crit == TIME)   return 1.0 + std::max(0.0, P.time_p - thresh) * gamma;
				if (crit == VOLUME) return 1.0 + std::max(0.0, P.volume_p - thresh) * gamma;
				if (crit == WEIGHT) return 1.0 + std::max(0.0, P.weight_p - thresh) * gamma;
				return 1.0; // SCORE neutral
			};

		auto perf = [&](int crit) -> double
			{
				return crit_stats[crit].ema_gain + EPS_GAIN;
			};

		std::array<double, N_CRIT> w =
		{
			perf(SCORE) * ctx_mult(SCORE),
			perf(TIME) * ctx_mult(TIME),
			perf(VOLUME) * ctx_mult(VOLUME),
			perf(WEIGHT) * ctx_mult(WEIGHT)
		};

		// sanitize
		for (double& wi : w)
		{
			if (!(wi > 0.0) || !std::isfinite(wi)) wi = 0.0;
		}

		// normalize
		double maxw = *std::max_element(w.begin(), w.end());
		if (maxw > 0.0)
		{
			for (double& wi : w) wi /= maxw;
		}

		// mix with floor
		for (double& wi : w) wi = 0.85 * wi + 0.15;

		auto pick_weighted = [&](const std::array<double, N_CRIT>& ww) -> int
			{
				double sum = 0.0;
				for (double x : ww) sum += x;

				if (!(sum > 0.0))
				{
					return std::uniform_int_distribution<int>(0, N_CRIT - 1)(engine);
				}

				double u = std::uniform_real_distribution<double>(0.0, sum)(engine);
				double acc = 0.0;
				for (int i = 0; i < N_CRIT; ++i)
				{
					acc += ww[i];
					if (u < acc) return i;
				}
				return N_CRIT - 1;
			};

		int crit_id = -1;
		if (forced_crit >= 0 && forced_crit < N_CRIT)
		{
			crit_id = forced_crit;   // fixed criterion for ablation
		}
		else
		{
			crit_id = pick_weighted(w); // adaptive behavior
		}

		auto kind = static_cast<RatioKind>(crit_id);

		// select enabled replacement neighborhood at random
		int pick = enabled_nb[std::uniform_int_distribution<int>(0, (int)enabled_nb.size() - 1)(engine)];

		switch (pick)
		{
		case 1:
		{
			auto nb = one_one_replace_gen_nb(s, tabulist, gb.score, kind);
			if (nb.size() == 0)
			{
				++nonb1;
				++no_feasible_moves_in_a_row;
			}
			else
			{
				no_feasible_moves_in_a_row = 0;
				executeMove(nb, s, tabulist);
				moved = true;
				tabulist.nextIteration();
			}
			break;
		}
		case 2:
		{
			auto nb = two_one_replace_gen_nb(s, tabulist, gb.score, kind);
			if (nb.size() == 0)
			{
				++nonb2;
				++no_feasible_moves_in_a_row;
			}
			else
			{
				no_feasible_moves_in_a_row = 0;
				executeMove(nb, s, tabulist);
				moved = true;
				tabulist.nextIteration();
			}
			break;
		}
		case 3:
		{
			auto nb = one_two_replace_gen_nb(s, tabulist, gb.score, kind);
			if (nb.size() == 0)
			{
				++nonb3;
				++no_feasible_moves_in_a_row;
			}
			else
			{
				no_feasible_moves_in_a_row = 0;
				executeMove(nb, s, tabulist);
				moved = true;
				tabulist.nextIteration();
			}
			break;
		}
		default:
			throw std::runtime_error("Invalid neighborhood selector value in Tabu::solve.");
		}

		
		if (no_feasible_moves_in_a_row > 5)
		{
			int idx = elites.pick_idx(s, ins->maxvertices, engine);
			if (idx >= 0)
			{
				std::unordered_set<uint64_t> A; ElitePool::fill_arcs(s, ins->maxvertices, A);
				std::unordered_set<uint64_t> B; ElitePool::fill_arcs(elites.get(idx), ins->maxvertices, B);
				double dist = ElitePool::arc_distance_frac(A, B);
				restart_dist_sum += dist;
				++restart_count;
				s = elites.get(idx);
				elites.mark_used(idx);
			}
			else
			{
				s = gb;
			}
			tabulist.clear();
			no_feasible_moves_in_a_row = 0;
			continue;
		}
		

		// --- RVND ---
		enum class LsOp { TwoOpt, Swap, Swap2, Relocate };
		std::array<LsOp, 4> ops = { LsOp::TwoOpt, LsOp::Swap, LsOp::Swap2, LsOp::Relocate };

		bool improved_any = true;
		while (improved_any)
		{
			improved_any = false;
			std::shuffle(ops.begin(), ops.end(), engine);

			for (auto op : ops)
			{
				bool improved = false;
				switch (op)
				{
				case LsOp::TwoOpt:   improved = two_opt_nb(s, 1);   break;
				case LsOp::Swap:     improved = swap_nb(s, 1);      break;
				case LsOp::Swap2:    improved = swap2_nb(s, 1);     break;
				case LsOp::Relocate: improved = relocate_nb(s, 1);  break;
				}
				if (improved)
				{
					improved_any = true;
					break;
				}
			}
		}

		// update criterion stats only in adaptive mode
		if (forced_crit < 0)
		{
			if (moved)
			{
				double gain = std::max(0, s.score - prev_score);
				crit_stats[crit_id].ema_gain =
					(1.0 - ema_rho) * crit_stats[crit_id].ema_gain + ema_rho * gain;
				crit_stats[crit_id].used++;
			}
			else
			{
				for (int c = 0; c < N_CRIT; ++c)
					crit_stats[c].ema_gain *= (1.0 - 0.05 * ema_rho);
			}
		}

		if (s.score > gb.score)
		{
			gb = s;
			noimpr = 0;
		}
		else
		{
			++noimpr;
		}

		if ((iter & 2) == 0) elites.consider(s, ins->maxvertices);
		elites.tick();
		++iter;
	}

	double cpuTime = std::chrono::duration<double>(clock::now() - t0).count();
	cout << gb << endl;
	gb.check();

	std::cout << "iter without replacement nb: "
		<< nonb1 << " <> " << nonb2 << " <> " << nonb3
		<< " total iterations: " << iter << std::endl;

	return Res(gb, cpuTime, bestknown);
}


Alns::Alns(Ins& ins, int iter_max_best, double alpha, int segment_len, int iter_per_segment, double sigma1, double sigma2, double sigma3, double temp_factor,double temp_min, double beta_frac, double lambda) : Moves(ins), iter_max_best(iter_max_best), alpha(alpha), segment_len(segment_len), iter_per_segment(iter_per_segment), sigma1(sigma1), sigma2(sigma2), sigma3(sigma3), temp_factor(temp_factor), temp_min(temp_min), beta_frac(beta_frac), lambda(lambda)
{
	gb = Sol(ins);//best sol
	s = Sol(ins);//iter sol
}


std::vector<std::pair<int, int>> Alns::collect_removable_positions(Sol& sol)
{
	std::vector<std::pair<int, int>> pos;
	const int startDepot = 0;
	const int endDepot = ins->maxvertices - 1;

	for (int t = 0; t < ins->maxtours; ++t)
	{
		const Tour& tour = sol.tours[t];
		for (int p = 1; p < (int)tour.seq.size() - 1; ++p) // skip both depots
		{
			int vid = tour.seq[p]->index;
			if (vid != startDepot && vid != endDepot)
			{
				pos.push_back({ t, p });
			}
		}
	}
	return pos;
}

bool Alns::time_feasible(const Sol::Tour& tour) const
{
	constexpr double tol = 1e-6;

	const int n = static_cast<int>(tour.seq.size());
	if (n < 2 ||
		static_cast<int>(tour.action.size()) != n ||
		static_cast<int>(tour.deptime.size()) != n ||
		static_cast<int>(tour.max_shift.size()) != n)
		return false;

	const int r = tour.index;
	const auto& td = ins->t[r];
	const double EDT = td.EDT;
	const double B = ins->breakdur;

	if (tour.seq.front() != &ins->v[0] ||
		tour.seq.back() != &ins->v[ins->maxvertices - 1])
		return false;

	int breakCount = 0;
	int foundBreak = -1;
	for (int i = 0; i < n; ++i)
	{
		if (tour.action[i] == 1)
		{
			++breakCount;
			foundBreak = i;
		}
	}

	if (breakCount != 1 || foundBreak != tour.breakindex)
		return false;

	// Rebuild the schedule using the same break and time-window semantics
	// as Tour::check(), but return false silently on failure.
	double currenttime = EDT + tour.deptime[0];

	for (int i = 0; i < n - 1; ++i)
	{
		Ins::Vertex* next = tour.seq[i + 1];
		const double arrival = ins->arrival_time(
			tour.seq[i]->con[next->index], currenttime);

		double startService;

		if (tour.action[i + 1] == 1)
		{
			const bool endDepot =
				next->index == ins->maxvertices - 1;

			const double startBreak =
				(!endDepot && arrival < ins->breakstart)
				? ins->breakstart
				: arrival;

			if (startBreak > ins->breakend + tol)
				return false;

			startService = std::max(
				startBreak + B, next->LTW[r]);
		}
		else
		{
			startService = std::max(arrival, next->LTW[r]);
		}

		if (startService > next->UTW[r] + tol)
			return false;

		currenttime = startService + next->serv;

		// Ensure the stored forward schedule matches the rebuilt schedule.
		if (std::abs((currenttime - EDT) - tour.deptime[i + 1]) > tol)
			return false;
	}

	const double duration = currenttime - EDT;
	if (duration > td.T_max + tol)
		return false;

	// The insertion evaluator depends on these cached slack values.
	// Recompute them on a copy and reject stale or negative slack.
	Sol::Tour slackCheck = tour;
	slackCheck.calc_maxshift();

	for (int i = 1; i < n; ++i)
	{
		if (slackCheck.max_shift[i] < -tol)
			return false;

		if (std::abs(slackCheck.max_shift[i] - tour.max_shift[i]) > 1e-2)
			return false;
	}

	return true;
}

bool Alns::try_remove_positions(Sol& sol, int tourSlot, std::vector<int> positions)
{
	if (positions.empty() || tourSlot < 0 || tourSlot >= ins->maxtours)
		return false;

	std::sort(positions.begin(), positions.end());
	positions.erase(
		std::unique(positions.begin(), positions.end()), positions.end());

	const int n = static_cast<int>(sol.tours[tourSlot].seq.size());
	for (int p : positions)
	{
		if (p <= 0 || p >= n - 1)
			return false; // never remove either depot
	}

	Sol trial = sol;

	if (positions.size() == 1)
	{
		trial.remove_vertex(trial.tours[tourSlot], positions[0]);
	}
	else if (positions.size() == 2)
	{
		trial.remove_vertices(
			trial.tours[tourSlot], positions[0], positions[1]);
	}
	else
	{
		// Delete from back to front so the remaining indices stay valid.
		for (auto it = positions.rbegin(); it != positions.rend(); ++it)
			trial.remove_vertex(trial.tours[tourSlot], *it);
	}

	if (!time_feasible(trial.tours[tourSlot]))
		return false;

	sol = std::move(trial);
	return true;
}

bool Alns::try_remove_vertex(Sol& sol, int tourSlot, int position)
{
	return try_remove_positions(sol, tourSlot, { position });
}

bool Alns::try_ranked_removal(
	Sol& sol,
	std::vector<RemovalCandidate> candidates,
	bool preferLargest,
	std::vector<RemovedCustomer>& removed)
{
	std::stable_sort(candidates.begin(), candidates.end(),
		[preferLargest](const RemovalCandidate& a,
			const RemovalCandidate& b)
		{
			return preferLargest
				? a.priority > b.priority
				: a.priority < b.priority;
		});

	// Failed attempts leave sol unchanged, so candidate positions remain valid.
	for (const RemovalCandidate& c : candidates)
	{
		const Tour& tour = sol.tours[c.tourSlot];
		if (c.position <= 0 ||
			c.position >= static_cast<int>(tour.seq.size()) - 1)
			continue;

		Ins::Vertex* v = tour.seq[c.position];
		const int pred = tour.seq[c.position - 1]->index;
		const int succ = tour.seq[c.position + 1]->index;

		if (try_remove_vertex(sol, c.tourSlot, c.position))
		{
			removed.push_back(
				{ v, c.tourSlot, c.position, pred, succ });
			return true;
		}
	}

	return false;
}

std::vector<RemovedCustomer> Alns::random_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	while (static_cast<int>(removed.size()) < beta)
	{
		auto positions = collect_removable_positions(sol);
		std::shuffle(positions.begin(), positions.end(), engine);

		bool removedOne = false;

		for (const auto& [tourSlot, position] : positions)
		{
			const Tour& tour = sol.tours[tourSlot];
			Ins::Vertex* v = tour.seq[position];
			const int pred = tour.seq[position - 1]->index;
			const int succ = tour.seq[position + 1]->index;

			if (try_remove_vertex(sol, tourSlot, position))
			{
				removed.push_back(
					{ v, tourSlot, position, pred, succ });
				removedOne = true;
				break;
			}
		}

		if (!removedOne)
			break;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::worst_remove_burden(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	constexpr double lambda_t = 2.0;
	constexpr double lambda_w = 1.0;
	constexpr double lambda_v = 1.0;

	for (int k = 0; k < beta; ++k)
	{
		std::vector<RemovalCandidate> candidates;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];
			const int route = tour.index;
			const double Wmax = ins->t[route].W_max;
			const double Vmax = ins->t[route].V_max;
			const double Tmax = ins->t[route].T_max;

			for (int p = 1; p < static_cast<int>(tour.seq.size()) - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];
				Ins::Vertex* x = tour.seq[p - 1];
				Ins::Vertex* z = tour.seq[p + 1];

				double detour =
					x->con[y->index]->determin +
					y->serv +
					y->con[z->index]->determin -
					x->con[z->index]->determin;

				detour = std::max(0.0, detour);

				const double timeFrac =
					Tmax > 0.0 ? detour / Tmax : 0.0;
				const double weightFrac =
					Wmax > 0.0 ? y->weight / Wmax : 0.0;
				const double volumeFrac =
					Vmax > 0.0 ? y->volume / Vmax : 0.0;

				const double burden =
					lambda_t * timeFrac +
					lambda_w * weightFrac +
					lambda_v * volumeFrac;

				const double value =
					double(y->score) / std::max(1e-9, burden);

				candidates.push_back({ value, t, p });
			}
		}

		// The original operator removes the smallest score-to-burden value.
		if (!try_ranked_removal(sol, std::move(candidates), false, removed))
			break;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::largest_time_saving_remove(Sol& sol, int beta)
{
	struct Candidate
	{
		double saving;
		int tourSlot;
		int position;
	};

	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		std::vector<Candidate> candidates;

		// Rank the current customers by estimated time saving.
		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];

			for (int p = 1; p < static_cast<int>(tour.seq.size()) - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];
				Ins::Vertex* x = tour.seq[p - 1];
				Ins::Vertex* z = tour.seq[p + 1];

				double saving =
					x->con[y->index]->determin
					+ y->serv
					+ y->con[z->index]->determin
					- x->con[z->index]->determin;

				candidates.push_back({ saving, t, p });
			}
		}

		std::sort(candidates.begin(), candidates.end(),
			[](const Candidate& a, const Candidate& b)
			{
				return a.saving > b.saving;
			});

		bool removed_one = false;

		// Try candidates from largest estimated saving to smallest.
		for (const Candidate& c : candidates)
		{
			// Failed attempts leave sol unchanged, so these positions remain valid.
			const Tour& tour = sol.tours[c.tourSlot];
			Ins::Vertex* v = tour.seq[c.position];
			const int pred = tour.seq[c.position - 1]->index;
			const int succ = tour.seq[c.position + 1]->index;

			if (try_remove_vertex(sol, c.tourSlot, c.position))
			{
				removed.push_back(
					{ v, c.tourSlot, c.position, pred, succ });
				removed_one = true;
				break;
			}
		}

		// No single-customer removal produced a feasible solution.
		if (!removed_one)
			break;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::sequence_remove(Sol& sol, int beta)
{
	struct Block
	{
		int tourSlot;
		int start;
		int length;
	};

	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	while (static_cast<int>(removed.size()) < beta)
	{
		const int remaining = beta - static_cast<int>(removed.size());
		std::vector<Block> blocks;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const int nCustomers =
				static_cast<int>(sol.tours[t].seq.size()) - 2;
			const int maxLen = std::min({ 3, remaining, nCustomers });

			for (int len = 1; len <= maxLen; ++len)
			{
				const int lastStart =
					static_cast<int>(sol.tours[t].seq.size()) - 1 - len;

				for (int start = 1; start <= lastStart; ++start)
					blocks.push_back({ t, start, len });
			}
		}

		std::shuffle(blocks.begin(), blocks.end(), engine);

		bool removedBlock = false;

		for (const Block& block : blocks)
		{
			const Tour& tour = sol.tours[block.tourSlot];
			std::vector<int> positions;
			std::vector<RemovedCustomer> blockInfo;

			for (int p = block.start; p < block.start + block.length; ++p)
			{
				positions.push_back(p);
				blockInfo.push_back({
					tour.seq[p],
					block.tourSlot,
					p,
					tour.seq[p - 1]->index,
					tour.seq[p + 1]->index
					});
			}

			if (try_remove_positions(sol, block.tourSlot, positions))
			{
				removed.insert(
					removed.end(), blockInfo.begin(), blockInfo.end());
				removedBlock = true;
				break;
			}
		}

		if (!removedBlock)
			break;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::break_neighborhood_remove(
	Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	std::vector<int> candidateTours;
	for (int t = 0; t < ins->maxtours; ++t)
		if (sol.tours[t].seq.size() > 2)
			candidateTours.push_back(t);

	if (candidateTours.empty()) return removed;

	std::uniform_int_distribution<int> distTour(
		0, static_cast<int>(candidateTours.size()) - 1);
	const int t = candidateTours[distTour(engine)];
	const Tour& tour = sol.tours[t];
	const int n = static_cast<int>(tour.seq.size());

	int b = tour.breakindex;
	if (b < 1 || b >= n) b = n - 1;

	// Order positions by proximity to the break.
	std::vector<int> ordered;
	if (b >= 1 && b <= n - 2)
		ordered.push_back(b);

	for (int radius = 1; ordered.size() < static_cast<size_t>(n - 2); ++radius)
	{
		const int left = b - radius;
		const int right = b + radius;

		if (left >= 1) ordered.push_back(left);
		if (right <= n - 2) ordered.push_back(right);

		if (left < 1 && right > n - 2)
			break;
	}

	const int desired = std::min(beta, static_cast<int>(ordered.size()));

	for (int len = desired; len >= 1; --len)
	{
		std::vector<int> positions(
			ordered.begin(), ordered.begin() + len);
		std::vector<RemovedCustomer> blockInfo;

		for (int p : positions)
		{
			blockInfo.push_back({
				tour.seq[p],
				t,
				p,
				tour.seq[p - 1]->index,
				tour.seq[p + 1]->index
				});
		}

		if (try_remove_positions(sol, t, positions))
		{
			removed = std::move(blockInfo);
			return removed;
		}
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::largest_demand_remove(
	Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		std::vector<RemovalCandidate> candidates;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];
			const int route = tour.index;
			const double Wmax = ins->t[route].W_max;
			const double Vmax = ins->t[route].V_max;

			for (int p = 1; p < static_cast<int>(tour.seq.size()) - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];
				const double demand =
					(Wmax > 0.0 ? y->weight / Wmax : 0.0) +
					(Vmax > 0.0 ? y->volume / Vmax : 0.0);

				candidates.push_back({ demand, t, p });
			}
		}

		if (!try_ranked_removal(sol, std::move(candidates), true, removed))
			break;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::lowest_profit_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		std::vector<RemovalCandidate> candidates;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];

			for (int p = 1; p < static_cast<int>(tour.seq.size()) - 1; ++p)
				candidates.push_back(
					{ double(tour.seq[p]->score), t, p });
		}

		if (!try_ranked_removal(sol, std::move(candidates), false, removed))
			break;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::largest_service_time_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		std::vector<RemovalCandidate> candidates;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];

			for (int p = 1; p < static_cast<int>(tour.seq.size()) - 1; ++p)
				candidates.push_back(
					{ tour.seq[p]->serv, t, p });
		}

		if (!try_ranked_removal(sol, std::move(candidates), true, removed))
			break;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::random_route_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;

	std::vector<int> candidateTours;
	for (int t = 0; t < ins->maxtours; ++t)
	{
		if ((int)sol.tours[t].seq.size() > 2) // has at least one customer
			candidateTours.push_back(t);
	}

	if (candidateTours.empty()) return removed;

	std::uniform_int_distribution<int> distTour(0, (int)candidateTours.size() - 1);
	int t = candidateTours[distTour(engine)];
	Tour& tour = sol.tours[t];

	// collect metadata first
	for (int p = 1; p < (int)tour.seq.size() - 1; ++p)
	{
		Ins::Vertex* v = tour.seq[p];
		int pred = tour.seq[p - 1]->index;
		int succ = tour.seq[p + 1]->index;
		removed.push_back({ v, t, p, pred, succ });
	}

	// remove all customers from back to front
	for (int p = (int)tour.seq.size() - 2; p >= 1; --p)
	{
		sol.remove_vertex(sol.tours[t], p);
	}

	return removed;
}

std::vector<Ins::Vertex*> Alns::collect_available_customers(const Sol& sol) const
{
	std::vector<Ins::Vertex*> pool;
	for (int vid = 1; vid < ins->maxvertices - 1; ++vid)
	{
		if (sol.available[vid])
			pool.push_back(&ins->v[vid]);
	}
	return pool;
}

Ins::Vertex* Alns::random_selection_from_pool(const std::vector<Ins::Vertex*>& pool) const
{
	if (pool.empty()) return nullptr;
	std::uniform_int_distribution<int> dist(0, (int)pool.size() - 1);
	return pool[dist(engine)];
}

Ins::Vertex* Alns::highest_score_selection_from_pool(const std::vector<Ins::Vertex*>& pool) const
{
	if (pool.empty()) return nullptr;

	return *std::max_element(pool.begin(), pool.end(),
		[](Ins::Vertex* a, Ins::Vertex* b)
		{
			return a->score < b->score;
		});
}

Ins::Vertex* Alns::highest_score_to_burden_selection_from_pool(	const std::vector<Ins::Vertex*>& pool) const
{
	if (pool.empty()) return nullptr;

	Ins::Vertex* best = nullptr;
	double bestVal = -DBL_MAX;

	double avgW = 0.0, avgV = 0.0;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		avgW += ins->t[d].W_max;
		avgV += ins->t[d].V_max;
	}
	avgW /= std::max(1, ins->maxtours);
	avgV /= std::max(1, ins->maxtours);

	for (auto* y : pool)
	{
		double burden =
			(avgW > 0.0 ? y->weight / avgW : 0.0) +
			(avgV > 0.0 ? y->volume / avgV : 0.0) +
			0.01 * y->serv;

		double val = double(y->score) / std::max(1e-9, burden);

		if (val > bestVal)
		{
			bestVal = val;
			best = y;
		}
	}
	return best;
}

Ins::Vertex* Alns::lowest_resource_selection_from_pool(	const std::vector<Ins::Vertex*>& pool) const
{
	if (pool.empty()) return nullptr;

	double avgW = 0.0, avgV = 0.0;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		avgW += ins->t[d].W_max;
		avgV += ins->t[d].V_max;
	}
	avgW /= std::max(1, ins->maxtours);
	avgV /= std::max(1, ins->maxtours);

	return *std::min_element(pool.begin(), pool.end(),
		[avgW, avgV](Ins::Vertex* a, Ins::Vertex* b)
		{
			double ra =
				(avgW > 0.0 ? a->weight / avgW : 0.0) +
				(avgV > 0.0 ? a->volume / avgV : 0.0);
			double rb =
				(avgW > 0.0 ? b->weight / avgW : 0.0) +
				(avgV > 0.0 ? b->volume / avgV : 0.0);
			return ra < rb;
		});
}


Ins::Vertex* Alns::dynamic_travel_time_profit_selection_from_pool(const Sol& sol,const std::vector<Ins::Vertex*>& pool,	const std::vector<RemovedCustomer>& removed) const
{
	if (pool.empty()) return nullptr;

	Ins::Vertex* bestNode = nullptr;
	double bestPref = -DBL_MAX;

	double maxScore = 1.0;
	for (auto* y : pool)
		maxScore = std::max(maxScore, double(y->score));

	double M = 1.0 + std::max(1.0, ins->t[0].T_max);
	for (int d = 1; d < ins->maxtours; ++d)
		M = std::max(M, 1.0 + ins->t[d].T_max);

	for (auto* y : pool)
	{
		std::vector<double> shifts = collect_feasible_insertion_shifts(sol, y, removed);
		if (shifts.empty()) continue;

		double u = std::uniform_real_distribution<double>(0.0, 1.0)(engine);
		if (u <= 1e-12) u = 1e-12;

		double profitScale = double(y->score) / maxScore;
		double pref = -DBL_MAX;

		if (shifts.size() == 1)
			pref = (M - shifts[0]) * u * profitScale;
		else
			pref = (shifts[1] - shifts[0]) * u * profitScale;

		if (pref > bestPref)
		{
			bestPref = pref;
			bestNode = y;
		}
	}

	return bestNode;
}

std::vector<double> Alns::collect_feasible_insertion_shifts(
	const Sol& sol,
	Ins::Vertex* y,
	const std::vector<RemovedCustomer>& removed) const
{
	std::vector<double> shifts;

	for (int d = 0; d < ins->maxtours; ++d)
	{
		const Sol::Tour& tour = sol.tours[d];
		for (int j = 0; j < (int)tour.seq.size() - 1; ++j)
		{
			GreedyInsertion cand;
			if (evaluate_insertion_position(sol, y, d, j, removed, cand))
			{
				shifts.push_back(cand.shift);
			}
		}
	}

	std::sort(shifts.begin(), shifts.end());
	return shifts;
}


Ins::Vertex* Alns::apply_selection_operator_from_pool(const Sol& sol,const std::vector<Ins::Vertex*>& pool,	SelectionOp sel_op,	const std::vector<RemovedCustomer>& removed) const
{
	switch (sel_op)
	{
	case SelectionOp::Random:
		return random_selection_from_pool(pool);
	case SelectionOp::HighestScore:
		return highest_score_selection_from_pool(pool);
	case SelectionOp::HighestScoreToBurden:
		return highest_score_to_burden_selection_from_pool(pool);
	case SelectionOp::LowestResource:
		return lowest_resource_selection_from_pool(pool);
	case SelectionOp::RegretTime:
		return dynamic_travel_time_profit_selection_from_pool(sol, pool, removed);
	}
	return nullptr;
}

void Alns::repair_with_selection_and_insertion_improved(Sol& sol,SelectionOp sel_op,InsertionOp ins_op,	const std::vector<RemovedCustomer>& removed)
{
	while (true)
	{
		std::vector<Ins::Vertex*> pool = collect_available_customers(sol);
		if (pool.empty()) break;

		bool inserted_any = false;

		while (!pool.empty())
		{
			Ins::Vertex* y = apply_selection_operator_from_pool(sol, pool, sel_op, removed);
			if (y == nullptr) break;

			GreedyInsertion insCand = apply_insertion_operator(sol, y, ins_op, removed);

			if (insCand.feasible)
			{
				sol.insert_vertex(sol.tours[insCand.tour_idx], y, insCand.position);
				inserted_any = true;
				break; // restart with refreshed pool
			}

			// remove y from temporary pool and try another node
			pool.erase(
				std::remove_if(pool.begin(), pool.end(),
					[y](Ins::Vertex* v) { return v->index == y->index; }),
				pool.end());
		}

		if (!inserted_any) break;
	}
}


std::vector<Ins::Vertex*> Alns::collect_insertable_customers(const Sol& sol,InsertionOp ins_op,	const std::vector<RemovedCustomer>& removed) const
{
	std::vector<Ins::Vertex*> pool;

	for (int vid = 1; vid < ins->maxvertices - 1; ++vid)
	{
		if (!sol.available[vid]) continue;

		Ins::Vertex* y = &ins->v[vid];
		GreedyInsertion cand = apply_insertion_operator(sol, y, ins_op, removed);

		if (cand.feasible)
			pool.push_back(y);
	}

	return pool;
}


const RemovedCustomer* Alns::find_removed_info(	Ins::Vertex* y,	const std::vector<RemovedCustomer>& removed) const
{
	for (const auto& r : removed)
	{
		if (r.v->index == y->index)
			return &r;
	}
	return nullptr;
}

void Alns::repair_with_selection_and_insertion(	Sol& sol,	SelectionOp sel_op,	InsertionOp ins_op,	const std::vector<RemovedCustomer>& removed)
{
	while (true)
	{
		// only nodes that can actually be inserted with the chosen insertion operator
		std::vector<Ins::Vertex*> pool =
			collect_insertable_customers(sol, ins_op, removed);

		if (pool.empty())
			break;

		// choose one insertable node according to the selected strategy
		Ins::Vertex* y = apply_selection_operator_from_pool(sol, pool, sel_op, removed);
		if (y == nullptr)
			break;

		// place it with the chosen insertion operator
		GreedyInsertion insCand = apply_insertion_operator(sol, y, ins_op, removed);
		if (!insCand.feasible)
			break; // should normally not happen because pool was prefiltered

		sol.insert_vertex(sol.tours[insCand.tour_idx], y, insCand.position);
	}
}

GreedyInsertion Alns::best_position_insertion(const Sol& sol,Ins::Vertex* y,const std::vector<RemovedCustomer>& removed) const
{
	GreedyInsertion best;

	for (int d = 0; d < ins->maxtours; ++d)
	{
		const Sol::Tour& tour = sol.tours[d];
		for (int j = 0; j < (int)tour.seq.size() - 1; ++j)
		{
			GreedyInsertion cand;
			if (!evaluate_insertion_position(sol, y, d, j, removed, cand)) continue;

			if (!best.feasible || cand.key > best.key)
				best = cand;
		}
	}
	return best;
}

GreedyInsertion Alns::first_feasible_insertion(	const Sol& sol,	Ins::Vertex* y,	const std::vector<RemovedCustomer>& removed) const
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		const Sol::Tour& tour = sol.tours[d];
		for (int j = 0; j < (int)tour.seq.size() - 1; ++j)
		{
			GreedyInsertion cand;
			if (evaluate_insertion_position(sol, y, d, j, removed, cand))
				return cand;
		}
	}
	return GreedyInsertion{};
}

GreedyInsertion Alns::last_feasible_insertion(	const Sol& sol,	Ins::Vertex* y,	const std::vector<RemovedCustomer>& removed) const
{
	GreedyInsertion last;

	for (int d = 0; d < ins->maxtours; ++d)
	{
		const Sol::Tour& tour = sol.tours[d];
		for (int j = 0; j < (int)tour.seq.size() - 1; ++j)
		{
			GreedyInsertion cand;
			if (evaluate_insertion_position(sol, y, d, j, removed, cand))
				last = cand;
		}
	}
	return last;
}

GreedyInsertion Alns::random_feasible_insertion(	const Sol& sol,	Ins::Vertex* y,	const std::vector<RemovedCustomer>& removed) const
{
	std::vector<GreedyInsertion> feasible;

	for (int d = 0; d < ins->maxtours; ++d)
	{
		const Sol::Tour& tour = sol.tours[d];
		for (int j = 0; j < (int)tour.seq.size() - 1; ++j)
		{
			GreedyInsertion cand;
			if (evaluate_insertion_position(sol, y, d, j, removed, cand))
				feasible.push_back(cand);
		}
	}

	if (feasible.empty()) return GreedyInsertion{};

	std::uniform_int_distribution<int> dist(0, (int)feasible.size() - 1);
	return feasible[dist(engine)];
}

GreedyInsertion Alns::least_loaded_best_position_insertion(	const Sol& sol,	Ins::Vertex* y,	const std::vector<RemovedCustomer>& removed) const
{
	int bestTour = -1;
	double bestLoad = std::numeric_limits<double>::infinity();

	for (int d = 0; d < ins->maxtours; ++d)
	{
		const Sol::Tour& tour = sol.tours[d];
		const double Wmax = ins->t[d].W_max;
		const double Vmax = ins->t[d].V_max;

		if (tour.weight + y->weight > Wmax + 1e-9) continue;
		if (tour.volume + y->volume > Vmax + 1e-9) continue;

		double load =
			(Wmax > 0.0 ? tour.weight / Wmax : 0.0) +
			(Vmax > 0.0 ? tour.volume / Vmax : 0.0);

		if (load < bestLoad)
		{
			bestLoad = load;
			bestTour = d;
		}
	}

	if (bestTour < 0) return GreedyInsertion{};

	GreedyInsertion best;
	const Sol::Tour& tour = sol.tours[bestTour];
	for (int j = 0; j < (int)tour.seq.size() - 1; ++j)
	{
		GreedyInsertion cand;
		if (!evaluate_insertion_position(sol, y, bestTour, j, removed, cand)) continue;

		if (!best.feasible || cand.key > best.key)
			best = cand;
	}
	return best;
}

bool Alns::evaluate_insertion_position(
	const Sol& sol,
	Ins::Vertex* y,
	int d,
	int j,
	const std::vector<RemovedCustomer>& removed,
	GreedyInsertion& out) const
{
	const Sol::Tour& tour = sol.tours[d];
	const int t = tour.index;
	const auto& Td = ins->t[t];

	const double EDT = Td.EDT;
	const double Wmax = Td.W_max;
	const double Vmax = Td.V_max;
	const double breakdur = ins->breakdur;

	if (tour.weight + y->weight > Wmax + 1e-9) return false;
	if (tour.volume + y->volume > Vmax + 1e-9) return false;

	Ins::Vertex* x = tour.seq[j];
	Ins::Vertex* z = tour.seq[j + 1];
	int breakz = tour.action[j + 1];

	if (!(x->nbi[t][y->index] && y->nbi[t][z->index])) return false;

	const std::vector<double> wait_suffix = tour.compute_wait_suffix();
	double currenttime = tour.deptime[j] + EDT;

	double base_up = ins->arrival_time(x->con[z->index], currenttime) - currenttime;
	double ins_lb = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
	double dt_lb = std::max(0.0, ins_lb - base_up);

	if (dt_lb > wait_suffix[j + 1] + tour.max_shift[j + 1] + 1e-9) return false;

	double at = ins->arrival_time(x->con[y->index], currenttime);

	if (at < y->LTW[t]) at = y->LTW[t];
	if (at > y->UTW[t] + 1e-9) return false;

	at += y->serv;
	at = ins->arrival_time(y->con[z->index], at);

	if (breakz)
	{
		const bool endp = (z->index == ins->maxvertices - 1);
		if (!endp && at < ins->breakstart) at = ins->breakstart;
		if (at > ins->breakend + 1e-9) return false;
		at += breakdur;
	}

	if (at < z->LTW[t]) at = z->LTW[t];
	at += z->serv;

	double shift = (at - EDT) - tour.deptime[j + 1];
	if (shift > tour.max_shift[j + 1] + 1e-9) return false;

	double rho =
		(Wmax > 0.0 ? y->weight / Wmax : 0.0) +
		(Vmax > 0.0 ? y->volume / Vmax : 0.0);

	double key = double(y->score) / std::max(1e-9, shift + rho);

	// penalize similar reinserts if recently removed
	const RemovedCustomer* rem = find_removed_info(y, removed);
	if (rem != nullptr)
	{
		double penalty = 1.0;
		int newPos = j + 1;

		if (d == rem->oldTour) penalty *= 0.90;
		if (newPos == rem->oldPos) penalty *= 0.75;
		if (x->index == rem->oldPred) penalty *= 0.65;
		if (z->index == rem->oldSucc) penalty *= 0.65;

		if (d == rem->oldTour &&
			x->index == rem->oldPred &&
			z->index == rem->oldSucc)
		{
			penalty *= 0.25;
		}
		key *= penalty;
	}

	out.feasible = true;
	out.vertex_idx = y->index;
	out.tour_idx = d;
	out.position = j;
	out.key = key;
	out.shift = shift;
	return true;

}

GreedyInsertion Alns::apply_insertion_operator(	const Sol& sol,	Ins::Vertex* y,	InsertionOp op,	const std::vector<RemovedCustomer>& removed) const
{
	switch (op)
	{
	case InsertionOp::BestPosition:
		return best_position_insertion(sol, y, removed);
	case InsertionOp::FirstFeasible:
		return first_feasible_insertion(sol, y, removed);
	case InsertionOp::LastFeasible:
		return last_feasible_insertion(sol, y, removed);
	case InsertionOp::RandomFeasible:
		return random_feasible_insertion(sol, y, removed);
	case InsertionOp::LeastLoadedBestPosition:
		return least_loaded_best_position_insertion(sol, y, removed);
	}
	return GreedyInsertion{};
}

bool Alns::accept_candidate(const Sol& cur, const Sol& cand, double T)
{
	if (cand.score >= cur.score) return true;

	double delta = double(cand.score - cur.score); // negative if worse
	double prob = std::exp(delta / std::max(1e-9, T));
	double u = std::uniform_real_distribution<double>(0.0, 1.0)(engine);

	return u < prob;
}

GreedyInsertion Alns::minimum_shift_insertion(
	const Sol& sol,
	Ins::Vertex* y,
	const std::vector<RemovedCustomer>& removed) const
{
	GreedyInsertion best;

	for (int d = 0; d < ins->maxtours; ++d)
	{
		const Sol::Tour& tour = sol.tours[d];

		for (int j = 0; j < static_cast<int>(tour.seq.size()) - 1; ++j)
		{
			GreedyInsertion candidate;
			if (!evaluate_insertion_position(sol, y, d, j, removed, candidate))
				continue;

			// Hammami describes choosing the best overall-time position.
			// Here, shift is the time-increase measure supplied by your evaluator.
			if (!best.feasible || candidate.shift < best.shift)
				best = candidate;
		}
	}

	return best;
}

bool Alns::ls_remove_one_and_refill(Sol& sol)
{
	auto positions = collect_removable_positions(sol);
	std::shuffle(positions.begin(), positions.end(), engine);

	for (const auto& [d, p] : positions)
	{
		const Tour& originalTour = sol.tours[d];
		Ins::Vertex* y = originalTour.seq[p];

		RemovedCustomer info{
			y,
			d,
			p,
			originalTour.seq[p - 1]->index,
			originalTour.seq[p + 1]->index
		};

		Sol candidate = sol;

		// The helper mutates candidate only if the resulting route is feasible.
		if (!try_remove_vertex(candidate, d, p))
			continue;

		const std::vector<RemovedCustomer> removed{ info };
		std::vector<int> touchedRoutes{ d };

		while (true)
		{
			auto pool = collect_available_customers(candidate);
			if (pool.empty())
				break;

			Ins::Vertex* next =
				dynamic_travel_time_profit_selection_from_pool(
					candidate, pool, removed);

			if (next == nullptr)
				break;

			GreedyInsertion place =
				minimum_shift_insertion(candidate, next, removed);

			if (!place.feasible)
				break;

			Sol afterInsert = candidate;
			afterInsert.insert_vertex(
				afterInsert.tours[place.tour_idx],
				next,
				place.position);

			if (!time_feasible(afterInsert.tours[place.tour_idx]))
				break;

			candidate = std::move(afterInsert);
			touchedRoutes.push_back(place.tour_idx);
		}

		bool feasible = true;
		for (int routeSlot : touchedRoutes)
		{
			if (!time_feasible(candidate.tours[routeSlot]))
			{
				feasible = false;
				break;
			}
		}

		if (feasible && candidate.score > sol.score)
		{
			sol = std::move(candidate);
			return true;
		}
	}

	return false;
}

bool Alns::ls_replace_with_higher_score(Sol& sol)
{
	auto positions = collect_removable_positions(sol);
	if (positions.empty())
		return false;

	std::uniform_int_distribution<int> pick(0, static_cast<int>(positions.size()) - 1);
	const auto [d, p] = positions[pick(engine)];

	Ins::Vertex* removedVertex = sol.tours[d].seq[p];
	RemovedCustomer info{
		removedVertex,
		d,
		p,
		sol.tours[d].seq[p - 1]->index,
		sol.tours[d].seq[p + 1]->index
	};

	Sol candidate = sol;
	candidate.remove_vertex(candidate.tours[d], p);

	const std::vector<RemovedCustomer> removed{ info };
	std::vector<Ins::Vertex*> pool;

	for (Ins::Vertex* y : collect_available_customers(candidate))
	{
		if (y->score <= removedVertex->score)
			continue;

		if (minimum_shift_insertion(candidate, y, removed).feasible)
			pool.push_back(y);
	}

	if (pool.empty())
		return false;

	Ins::Vertex* replacement = random_selection_from_pool(pool);
	GreedyInsertion place =
		random_feasible_insertion(candidate, replacement, removed);

	if (!place.feasible)
		return false;

	candidate.insert_vertex(
		candidate.tours[place.tour_idx], replacement, place.position);

	if (candidate.score > sol.score)
	{
		sol = candidate;
		return true;
	}

	return false;
}

bool Alns::ls_remove_two_insert_one(Sol& sol)
{
	auto positions = collect_removable_positions(sol);
	if (positions.size() < 2)
		return false;

	std::shuffle(positions.begin(), positions.end(), engine);

	auto make_removed_info =
		[&sol](const std::pair<int, int>& pos)
		{
			const int d = pos.first;
			const int p = pos.second;
			Ins::Vertex* y = sol.tours[d].seq[p];

			return RemovedCustomer{
				y,
				d,
				p,
				sol.tours[d].seq[p - 1]->index,
				sol.tours[d].seq[p + 1]->index
			};
		};

	for (size_t a = 0; a < positions.size(); ++a)
	{
		for (size_t b = a + 1; b < positions.size(); ++b)
		{
			const auto first = positions[a];
			const auto second = positions[b];

			std::vector<RemovedCustomer> removed{
				make_removed_info(first),
				make_removed_info(second)
			};

			Sol candidate = sol;

			// If both are on the same route, remove the later position first.
			if (first.first == second.first)
			{
				const int high = std::max(first.second, second.second);
				const int low = std::min(first.second, second.second);

				candidate.remove_vertex(candidate.tours[first.first], high);
				candidate.remove_vertex(candidate.tours[first.first], low);
			}
			else
			{
				candidate.remove_vertex(
					candidate.tours[first.first], first.second);
				candidate.remove_vertex(
					candidate.tours[second.first], second.second);
			}

			// Check the completed two-customer removal before using it as
			// the starting point for insertion evaluation.
			if (!time_feasible(candidate.tours[first.first]))
				continue;

			if (second.first != first.first &&
				!time_feasible(candidate.tours[second.first]))
				continue;

			std::vector<int> touchedRoutes{
				first.first, second.first
			};

			// Find an unserved customer other than the two just removed.
			auto pool = collect_available_customers(candidate);
			pool.erase(
				std::remove_if(
					pool.begin(), pool.end(),
					[&removed](Ins::Vertex* v)
					{
						return v->index == removed[0].v->index ||
							v->index == removed[1].v->index;
					}),
				pool.end());

			std::shuffle(pool.begin(), pool.end(), engine);

			bool insertedNewCustomer = false;

			for (Ins::Vertex* newCustomer : pool)
			{
				GreedyInsertion place =
					minimum_shift_insertion(candidate, newCustomer, removed);

				if (!place.feasible)
					continue;

				Sol afterInsert = candidate;
				afterInsert.insert_vertex(
					afterInsert.tours[place.tour_idx],
					newCustomer,
					place.position);

				if (!time_feasible(afterInsert.tours[place.tour_idx]))
					continue;

				candidate = std::move(afterInsert);
				touchedRoutes.push_back(place.tour_idx);
				insertedNewCustomer = true;
				break;
			}

			if (!insertedNewCustomer)
				continue;

			// Try to reinsert the two removed customers in random order.
			std::shuffle(removed.begin(), removed.end(), engine);

			for (const RemovedCustomer& oldCustomer : removed)
			{
				GreedyInsertion place =
					minimum_shift_insertion(
						candidate, oldCustomer.v, removed);

				if (!place.feasible)
					continue;

				Sol afterInsert = candidate;
				afterInsert.insert_vertex(
					afterInsert.tours[place.tour_idx],
					oldCustomer.v,
					place.position);

				if (!time_feasible(afterInsert.tours[place.tour_idx]))
					continue;

				candidate = std::move(afterInsert);
				touchedRoutes.push_back(place.tour_idx);
			}

			bool feasible = true;
			for (int routeSlot : touchedRoutes)
			{
				if (!time_feasible(candidate.tours[routeSlot]))
				{
					feasible = false;
					break;
				}
			}

			if (feasible && candidate.score > sol.score)
			{
				sol = std::move(candidate);
				return true;
			}
		}
	}

	return false;
}

void Alns::apply_local_search(Sol& sol)
{
	//two_opt_nb(sol, 1);
	//ls_remove_one_and_refill(sol);
	//ls_replace_with_higher_score(sol);
	//swap2_nb(sol, 1);
	//ls_remove_two_insert_one(sol);
}

void Alns::update_operator_weights(double lambda)
{
	for (int i = 0; i < 7; ++i)
	{
		if (rem_used[i] > 0)
			rem_w[i] = lambda * (rem_score[i] / rem_used[i]) + (1.0 - lambda) * rem_w[i];
		rem_score[i] = 0.0;
		rem_used[i] = 0;
	}

	for (int i = 0; i < 5; ++i)
	{
		if (sel_used[i] > 0)
			sel_w[i] = lambda * (sel_score[i] / sel_used[i]) + (1.0 - lambda) * sel_w[i];
		sel_score[i] = 0.0;
		sel_used[i] = 0;
	}

	for (int i = 0; i < 5; ++i)
	{
		if (ins_used[i] > 0)
			ins_w[i] = lambda * (ins_score[i] / ins_used[i]) + (1.0 - lambda) * ins_w[i];
		ins_score[i] = 0.0;
		ins_used[i] = 0;
	}
}

Res Alns::solve(int bestknown, double max_time_sec)
{
	using clock = std::chrono::steady_clock;
	auto t0 = clock::now();

	auto elapsed_sec = [&]() -> double
		{
			return std::chrono::duration<double>(clock::now() - t0).count();
		};

	auto time_up = [&]() -> bool
		{
			return elapsed_sec() >= max_time_sec;
		};

	s.reset();
	parallel_construct(s);
	gb = s;

	rem_w.fill(1.0); sel_w.fill(1.0); ins_w.fill(1.0);
	rem_score.fill(0.0); sel_score.fill(0.0); ins_score.fill(0.0);
	rem_used.fill(0); sel_used.fill(0); ins_used.fill(0);

	int noimpr = 0;
	double T = temp_factor;
	bool no_customers = false;

	for (int seg = 0;seg < segment_len && noimpr < iter_max_best && !no_customers;++seg)
	{
		//cout<<"seg: " << seg << endl;
		bool timed_out = false;

		for (int iter = 0;iter < iter_per_segment && noimpr < iter_max_best;++iter)
		{
			//cout<<"iter: " << iter << endl;
			if (time_up())
			{
				timed_out = true;
				//cout << "break in if time_up()" << endl;
				break;
			}

			// s is the last admissible solution; rejected candidates leave it unchanged.
			Sol cand = s;
			std::vector<RemovedCustomer> removed;

			int nVisited = 0;
			for (int d = 0; d < ins->maxtours; ++d)
				nVisited += std::max(
					0, static_cast<int>(cand.tours[d].seq.size()) - 2);

			if (nVisited == 0)
			{
				no_customers = true;
				//cout << "break due to no customers visited" << endl;
				break;
			}

			const int beta_max =
				std::max(1, static_cast<int>(std::floor(beta_frac * nVisited)));
			std::uniform_int_distribution<int> beta_picker(1, beta_max);
			const int beta = beta_picker(engine);

			// Select and apply a removal operator.
			const int rem_idx = weighted_pick(rem_w, engine);
			const int rem_op = rem_idx + 1;
			rem_used[rem_idx]++;

			switch (rem_op)
			{
			case 1: removed = random_remove(cand, beta); break;
			case 2: removed = largest_time_saving_remove(cand, beta); break;
			case 3: removed = largest_demand_remove(cand, beta); break;
			case 4: removed = lowest_profit_remove(cand, beta); break;
			case 5: removed = largest_service_time_remove(cand, beta); break;
			case 6: removed = random_route_remove(cand, beta); break;
			case 7: removed = sequence_remove(cand, beta); break;
			}
			cand.check();
			// Select node-selection and insertion operators independently.
			const int sel_idx = weighted_pick(sel_w, engine);
			const int ins_idx = weighted_pick(ins_w, engine);

			const SelectionOp sel_op =
				static_cast<SelectionOp>(sel_idx + 1);
			const InsertionOp ins_op =
				static_cast<InsertionOp>(ins_idx + 1);

			sel_used[sel_idx]++;
			ins_used[ins_idx]++;

			repair_with_selection_and_insertion(cand, sel_op, ins_op, removed);

			const bool improving_current = cand.score > s.score;
			const bool accepted = accept_candidate(s, cand, T);

			if (accepted && improving_current)
				apply_local_search(cand);

			const bool new_global_best =
				accepted && cand.score > gb.score;

			if (accepted)
				s = cand;

			if (new_global_best)
			{
				gb = cand;
				noimpr = 0;
			}
			else
			{
				++noimpr;
			}

			double reward = 0.0;
			if (new_global_best)
				reward = sigma1;
			else if (accepted && improving_current)
				reward = sigma2;
			else if (accepted)
				reward = sigma3;

			rem_score[rem_idx] += reward;
			sel_score[sel_idx] += reward;
			ins_score[ins_idx] += reward;

			// Hammami reheats at Tmin. The SPP call is intentionally omitted.
			if (T <= temp_min)
				T = temp_factor;

			T = std::max(temp_min, T * alpha);
		}

		// Update weights at the end of each completed segment.
		if (!timed_out && !no_customers)
			update_operator_weights(lambda);

		if (timed_out)
		{
			//cout << "break in bottom time_out" << endl;
			std::cerr << "timeout: elapsed=" << elapsed_sec()<< ", limit=" << max_time_sec << '\n';
			break;
		}
	}
	//cout <<"noimpr: " << noimpr << endl;
	double wallTime = elapsed_sec();
	//std::cerr << "solve return: wallTime=" << wallTime << '\n';
	gb.check();
	return Res(gb, wallTime, bestknown);
}

