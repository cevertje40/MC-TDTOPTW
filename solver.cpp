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


std::vector<RemovedCustomer> Alns::random_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; )
	{
		auto pos = collect_removable_positions(sol);
		if (pos.empty()) break;

		// if only one left or only one removal remaining
		if (pos.size() == 1 || beta - k == 1)
		{
			std::uniform_int_distribution<int> dist(0, (int)pos.size() - 1);
			auto [tour_idx, position] = pos[dist(engine)];

			Ins::Vertex* v = sol.tours[tour_idx].seq[position];
			int pred = sol.tours[tour_idx].seq[position - 1]->index;
			int succ = sol.tours[tour_idx].seq[position + 1]->index;

			removed.push_back({ v, tour_idx, position, pred, succ });
			sol.remove_vertex(sol.tours[tour_idx], position);
			++k;
			continue;
		}

		// try removing two at once if possible
		std::uniform_int_distribution<int> dist1(0, (int)pos.size() - 1);
		int idx1 = dist1(engine);

		std::uniform_int_distribution<int> dist2(0, (int)pos.size() - 2);
		int idx2 = dist2(engine);
		if (idx2 >= idx1) ++idx2;

		auto [tour1, pos1] = pos[idx1];
		auto [tour2, pos2] = pos[idx2];

		if (tour1 == tour2 && std::abs(pos1 - pos2) == 1)
		{
			if (pos1 > pos2) std::swap(pos1, pos2);

			Ins::Vertex* v1 = sol.tours[tour1].seq[pos1];
			Ins::Vertex* v2 = sol.tours[tour1].seq[pos2];

			int pred1 = sol.tours[tour1].seq[pos1 - 1]->index;
			int succ1 = sol.tours[tour1].seq[pos1 + 1]->index;

			int pred2 = sol.tours[tour1].seq[pos2 - 1]->index;
			int succ2 = sol.tours[tour1].seq[pos2 + 1]->index;

			removed.push_back({ v1, tour1, pos1, pred1, succ1 });
			removed.push_back({ v2, tour1, pos2, pred2, succ2 });

			sol.remove_vertices(sol.tours[tour1], pos1, pos2);
			k += 2;
		}
		else
		{
			Ins::Vertex* v = sol.tours[tour1].seq[pos1];
			int pred = sol.tours[tour1].seq[pos1 - 1]->index;
			int succ = sol.tours[tour1].seq[pos1 + 1]->index;

			removed.push_back({ v, tour1, pos1, pred, succ });
			sol.remove_vertex(sol.tours[tour1], pos1);
			++k;
		}
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::worst_remove_burden(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	const double lambda_t = 2.0;
	const double lambda_w = 1.0;
	const double lambda_v = 1.0;

	for (int k = 0; k < beta; ++k)
	{
		int bestTour = -1;
		int bestPos = -1;
		double worstValue = std::numeric_limits<double>::infinity();

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];
			const double Wmax = ins->t[t].W_max;
			const double Vmax = ins->t[t].V_max;
			const double Tmax = ins->t[t].T_max;

			for (int p = 1; p < (int)tour.seq.size() - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];
				Ins::Vertex* x = tour.seq[p - 1];
				Ins::Vertex* z = tour.seq[p + 1];

				double detour =
					x->con[y->index]->determin +
					y->serv +
					y->con[z->index]->determin -
					x->con[z->index]->determin;

				if (detour < 0.0) detour = 0.0;

				double timeFrac = (Tmax > 0.0 ? detour / Tmax : 0.0);
				double weightFrac = (Wmax > 0.0 ? y->weight / Wmax : 0.0);
				double volumeFrac = (Vmax > 0.0 ? y->volume / Vmax : 0.0);

				double burden = lambda_t * timeFrac
					+ lambda_w * weightFrac
					+ lambda_v * volumeFrac;

				double value = double(y->score) / std::max(1e-9, burden);

				if (value < worstValue)
				{
					worstValue = value;
					bestTour = t;
					bestPos = p;
				}
			}
		}

		if (bestTour < 0) break;

		Ins::Vertex* v = sol.tours[bestTour].seq[bestPos];
		int pred = sol.tours[bestTour].seq[bestPos - 1]->index;
		int succ = sol.tours[bestTour].seq[bestPos + 1]->index;

		removed.push_back({ v, bestTour, bestPos, pred, succ });
		sol.remove_vertex(sol.tours[bestTour], bestPos);
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::largest_time_saving_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		int bestTour = -1;
		int bestPos = -1;
		double bestSaving = -1.0;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];

			for (int p = 1; p < (int)tour.seq.size() - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];
				Ins::Vertex* x = tour.seq[p - 1];
				Ins::Vertex* z = tour.seq[p + 1];

				double saving =
					x->con[y->index]->determin +
					y->serv +
					y->con[z->index]->determin -
					x->con[z->index]->determin;

				if (saving > bestSaving)
				{
					bestSaving = saving;
					bestTour = t;
					bestPos = p;
				}
			}
		}

		if (bestTour < 0) break;

		Ins::Vertex* v = sol.tours[bestTour].seq[bestPos];
		int pred = sol.tours[bestTour].seq[bestPos - 1]->index;
		int succ = sol.tours[bestTour].seq[bestPos + 1]->index;

		removed.push_back({ v, bestTour, bestPos, pred, succ });
		sol.remove_vertex(sol.tours[bestTour], bestPos);
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::sequence_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	int removedCount = 0;

	while (removedCount < beta)
	{
		// routes with at least one customer
		std::vector<int> candidateTours;
		for (int t = 0; t < ins->maxtours; ++t)
		{
			if ((int)sol.tours[t].seq.size() > 2)
				candidateTours.push_back(t);
		}
		if (candidateTours.empty()) break;

		std::uniform_int_distribution<int> distTour(0, (int)candidateTours.size() - 1);
		int t = candidateTours[distTour(engine)];
		Tour& tour = sol.tours[t];

		int nCustomers = (int)tour.seq.size() - 2;
		if (nCustomers <= 0) break;

		int remaining = beta - removedCount;
		int maxLen = std::min(3, std::min(remaining, nCustomers));
		int minLen = std::min(2, maxLen);
		if (minLen <= 0) minLen = 1;

		std::uniform_int_distribution<int> distLen(minLen, maxLen);
		int len = distLen(engine);

		int lastStart = (int)tour.seq.size() - 1 - len;
		if (lastStart < 1) break;

		std::uniform_int_distribution<int> distStart(1, lastStart);
		int startPos = distStart(engine);

		// collect metadata before deletion
		for (int p = startPos; p < startPos + len; ++p)
		{
			Ins::Vertex* v = tour.seq[p];
			int pred = tour.seq[p - 1]->index;
			int succ = tour.seq[p + 1]->index;
			removed.push_back({ v, t, p, pred, succ });
		}

		// remove back to front
		for (int p = startPos + len - 1; p >= startPos; --p)
		{
			sol.remove_vertex(sol.tours[t], p);
		}

		removedCount += len;
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::break_neighborhood_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	// routes with customers
	std::vector<int> candidateTours;
	for (int t = 0; t < ins->maxtours; ++t)
	{
		if ((int)sol.tours[t].seq.size() > 2)
			candidateTours.push_back(t);
	}
	if (candidateTours.empty()) return removed;

	std::uniform_int_distribution<int> distTour(0, (int)candidateTours.size() - 1);
	int t = candidateTours[distTour(engine)];
	Tour& tour = sol.tours[t];

	int n = (int)tour.seq.size();
	if (n <= 2) return removed;

	int b = tour.breakindex;
	if (b < 1 || b >= n) b = n - 1;

	std::vector<int> candPos;

	// prioritize break position if internal
	if (b >= 1 && b <= n - 2) candPos.push_back(b);

	int left = b - 1;
	int right = b + 1;

	while ((int)candPos.size() < beta && (left >= 1 || right <= n - 2))
	{
		if (left >= 1)
		{
			candPos.push_back(left);
			--left;
			if ((int)candPos.size() >= beta) break;
		}
		if (right <= n - 2)
		{
			candPos.push_back(right);
			++right;
		}
	}

	if (candPos.empty())
	{
		// fallback: remove one random customer
		std::uniform_int_distribution<int> distPos(1, n - 2);
		int p = distPos(engine);

		Ins::Vertex* v = tour.seq[p];
		int pred = tour.seq[p - 1]->index;
		int succ = tour.seq[p + 1]->index;
		removed.push_back({ v, t, p, pred, succ });

		sol.remove_vertex(sol.tours[t], p);
		return removed;
	}

	std::sort(candPos.begin(), candPos.end());

	// collect metadata before deletion
	for (int p : candPos)
	{
		Ins::Vertex* v = tour.seq[p];
		int pred = tour.seq[p - 1]->index;
		int succ = tour.seq[p + 1]->index;
		removed.push_back({ v, t, p, pred, succ });
	}

	// remove from back to front
	for (int i = (int)candPos.size() - 1; i >= 0; --i)
	{
		sol.remove_vertex(sol.tours[t], candPos[i]);
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::largest_demand_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		int bestTour = -1;
		int bestPos = -1;
		double bestDemand = -1.0;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];
			const double Wmax = ins->t[t].W_max;
			const double Vmax = ins->t[t].V_max;

			for (int p = 1; p < (int)tour.seq.size() - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];

				double demand =
					(Wmax > 0.0 ? y->weight / Wmax : 0.0) +
					(Vmax > 0.0 ? y->volume / Vmax : 0.0);

				if (demand > bestDemand)
				{
					bestDemand = demand;
					bestTour = t;
					bestPos = p;
				}
			}
		}

		if (bestTour < 0) break;

		Ins::Vertex* v = sol.tours[bestTour].seq[bestPos];
		int pred = sol.tours[bestTour].seq[bestPos - 1]->index;
		int succ = sol.tours[bestTour].seq[bestPos + 1]->index;

		removed.push_back({ v, bestTour, bestPos, pred, succ });
		sol.remove_vertex(sol.tours[bestTour], bestPos);
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::lowest_profit_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		int bestTour = -1;
		int bestPos = -1;
		int lowestScore = std::numeric_limits<int>::max();

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];

			for (int p = 1; p < (int)tour.seq.size() - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];

				if (y->score < lowestScore)
				{
					lowestScore = y->score;
					bestTour = t;
					bestPos = p;
				}
			}
		}

		if (bestTour < 0) break;

		Ins::Vertex* v = sol.tours[bestTour].seq[bestPos];
		int pred = sol.tours[bestTour].seq[bestPos - 1]->index;
		int succ = sol.tours[bestTour].seq[bestPos + 1]->index;

		removed.push_back({ v, bestTour, bestPos, pred, succ });
		sol.remove_vertex(sol.tours[bestTour], bestPos);
	}

	return removed;
}

std::vector<RemovedCustomer> Alns::largest_service_time_remove(Sol& sol, int beta)
{
	std::vector<RemovedCustomer> removed;
	if (beta <= 0) return removed;

	for (int k = 0; k < beta; ++k)
	{
		int bestTour = -1;
		int bestPos = -1;
		double bestServ = -1.0;

		for (int t = 0; t < ins->maxtours; ++t)
		{
			const Tour& tour = sol.tours[t];

			for (int p = 1; p < (int)tour.seq.size() - 1; ++p)
			{
				Ins::Vertex* y = tour.seq[p];

				if (y->serv > bestServ)
				{
					bestServ = y->serv;
					bestTour = t;
					bestPos = p;
				}
			}
		}

		if (bestTour < 0) break;

		Ins::Vertex* v = sol.tours[bestTour].seq[bestPos];
		int pred = sol.tours[bestTour].seq[bestPos - 1]->index;
		int succ = sol.tours[bestTour].seq[bestPos + 1]->index;

		removed.push_back({ v, bestTour, bestPos, pred, succ });
		sol.remove_vertex(sol.tours[bestTour], bestPos);
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
	const auto& Td = ins->t[d];

	const double EDT = Td.EDT;
	const double Wmax = Td.W_max;
	const double Vmax = Td.V_max;
	const double breakdur = ins->breakdur;

	if (tour.weight + y->weight > Wmax + 1e-9) return false;
	if (tour.volume + y->volume > Vmax + 1e-9) return false;

	Ins::Vertex* x = tour.seq[j];
	Ins::Vertex* z = tour.seq[j + 1];
	int breakz = tour.action[j + 1];

	if (!(x->nbi[d][y->index] && y->nbi[d][z->index])) return false;

	const std::vector<double> wait_suffix = tour.compute_wait_suffix();
	double currenttime = tour.deptime[j] + EDT;

	double base_up = ins->arrival_time(x->con[z->index], currenttime) - currenttime;
	double ins_lb = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
	double dt_lb = std::max(0.0, ins_lb - base_up);

	if (dt_lb > wait_suffix[j + 1] + tour.max_shift[j + 1] + 1e-9) return false;

	double at = ins->arrival_time(x->con[y->index], currenttime);

	if (at < y->LTW[d]) at = y->LTW[d];
	if (at > y->UTW[d] + 1e-9) return false;

	at += y->serv;
	at = ins->arrival_time(y->con[z->index], at);

	if (breakz)
	{
		const bool endp = (z->index == ins->maxvertices - 1);
		if (!endp && at < ins->breakstart) at = ins->breakstart;
		if (at > ins->breakend + 1e-9) return false;
		at += breakdur;
	}

	if (at < z->LTW[d]) at = z->LTW[d];
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

void Alns::apply_local_search(Sol& sol)
{
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
			case LsOp::TwoOpt:   improved = two_opt_nb(sol, 1);   break;
			case LsOp::Swap:     improved = swap_nb(sol, 1);      break;
			case LsOp::Swap2:    improved = swap2_nb(sol, 1);     break;
			case LsOp::Relocate: improved = relocate_nb(sol, 1);  break;
			}
			if (improved)
			{
				improved_any = true;
				break;
			}
		}
	}
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

	int total_it = segment_len * iter_per_segment;
	// SA temperature
	double T = temp_factor;
	
	for (int iter = 0; iter < total_it; ++iter)
	{
		if (time_up())
		{
			break;
		}
		Sol cand = s;
		std::vector<RemovedCustomer> removed;

		// ----- destruction size beta -----
		int nVisited = 0;
		for (int t = 0; t < ins->maxtours; ++t)
			nVisited += std::max(0, (int)cand.tours[t].seq.size() - 2);

		if (nVisited == 0)
			break;

		// slightly more conservative than Hammami for your richer problem
		int beta_max = std::max(1, (int)std::floor(beta_frac * nVisited));
		std::uniform_int_distribution<int> beta_picker(1, beta_max);
		int beta = beta_picker(engine);

		// ----- choose removal operator -----

		int rem_idx = weighted_pick(rem_w, engine); 
		int rem_op = rem_idx + 1;
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

		// ----- choose selection + insertion operators -----
		
		int sel_idx = weighted_pick(sel_w, engine); // 0..4
		int ins_idx = weighted_pick(ins_w, engine); // 0..4

		SelectionOp sel_op = static_cast<SelectionOp>(sel_idx + 1);
		InsertionOp ins_op = static_cast<InsertionOp>(ins_idx + 1);

		
		sel_used[sel_idx]++;
		ins_used[ins_idx]++;

		// ----- repair -----
		repair_with_selection_and_insertion(cand, sel_op, ins_op, removed);

		// optional debug
		/*
		if (!cand.check())
		{
			std::cout << "Error after repair" << std::endl;
			std::abort();
		}
		*/

		// ----- acceptance -----
		bool improving_current = (cand.score > s.score);
		bool accepted = accept_candidate(s, cand, T);

		if (accepted && improving_current)
		{
			//apply_local_search(cand);
		}

		bool new_global_best = (cand.score > gb.score);

		if (accepted)
		{
			s = cand;
		}

		// ----- update global best -----
		if (new_global_best)
		{
			gb = cand;
			noimpr = 0;
		}
		else
		{
			++noimpr;
		}

		// ----- operator reward -----
		double reward = 0.0;
		if (new_global_best) reward = sigma1;
		else if (accepted && improving_current) reward = sigma2;
		else if (accepted) reward = sigma3;

		rem_score[rem_idx] += reward;
		sel_score[sel_idx] += reward;
		ins_score[ins_idx] += reward;

		// ----- cool temperature -----
		T = std::max(temp_min, T * alpha);

		// ----- simple restart / reheat -----
		if (noimpr > iter_max_best)
		{
			s = gb;
			noimpr = 0;
			T = temp_factor;
		}

		// ----- segment-based weight update -----
		if ((iter + 1) % iter_per_segment == 0)
		{
			update_operator_weights(lambda);
		}
	}

	double cpuTime = std::chrono::duration<double>(clock::now() - t0).count();
	gb.check();
	return Res(gb, cpuTime, bestknown);
}

