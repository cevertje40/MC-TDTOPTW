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


Tabu::Tabu(Ins& ins, int max_noimpr, int nb_tabu_it, double ema_rho, double thresh, double gamma) : Moves(ins), max_noimpr(max_noimpr), nb_tabu_it(nb_tabu_it), EMA_RHO(ema_rho), THRESH(thresh), GAMMA(gamma)
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

void Tabu::parallel_construct(Sol& sol)
{
	// 0) Per-tour depot-reachable lists, sorted by score
	std::vector<std::vector<Ins::Vertex*>> dep_reach(ins->maxtours);
	for (int d = 0; d < ins->maxtours; ++d) 
	{
		for (int i = 1; i < ins->maxvertices - 1; ++i) 
		{
			const double lb = ins->v[0].con[i]->determin+ ins->v[i].serv+ ins->v[i].con[ins->maxvertices - 1]->determin;
			if (lb <= ins->t[d].T_max + 1e-9)
				dep_reach[d].push_back(&ins->v[i]);
		}
		std::sort(dep_reach[d].begin(), dep_reach[d].end(),[](Ins::Vertex* a, Ins::Vertex* b) { return a->score > b->score; });
	}

	// Helper: try to append cand to the end of tour d; return new deptime, break flag
	auto try_append = [&](int d, Ins::Vertex* cand, double& new_deptime, int& brk) -> bool 
		{
			if (!sol.available[cand->index]) return false;

			Sol::Tour& tour = sol.tours[d];
			if (tour.weight + cand->weight > ins->t[d].W_max + 1e-9) return false;
			if (tour.volume + cand->volume > ins->t[d].V_max + 1e-9) return false;

			double currenttime = tour.deptime.back() + ins->t[d].EDT;
			Ins::Vertex* last = tour.seq.back();

			double arr = ins->arrival_time(last->con[cand->index], currenttime);
			int need_break = -1;
			if ((tour.breakindex == -1) &&((arr >= ins->breakstart) || (cand->LTW[d] - ins->breakdur >= ins->breakstart))) 
			{
				arr += ins->breakdur;
				need_break = 1;
			}
			else 
			{
				need_break = 0;
			}
			if (arr < cand->LTW[d]) arr = cand->LTW[d];
			if (arr > cand->UTW[d] + 1e-9) return false;

			arr += cand->serv;

			double enddep = ins->arrival_time(cand->con[ins->maxvertices - 1], arr);
			if ((tour.breakindex == -1) && (need_break == 0)) enddep += ins->breakdur;
			if (enddep > ins->t[d].LAT + 1e-9) return false;

			new_deptime = arr - ins->t[d].EDT; // includes break if taken
			brk = need_break;
			return true;
		};

	// 1) Phase A: seed each tour with one vertex (round-robin)
	std::vector<size_t> ptr(ins->maxtours, 0);   // per-tour scan pointer
	bool seeded_any = true;
	const int K = 12; // small scan window per pass to stay fast

	while (seeded_any) 
	{
		seeded_any = false;

		for (int d = 0; d < ins->maxtours; ++d) 
		{
			Sol::Tour& tour = sol.tours[d];

			double best_key = -DBL_MAX;
			Ins::Vertex* best_v = nullptr;
			double best_deptime = 0.0;
			int best_brk = 0;

			int examined = 0;
			auto& L = dep_reach[d];
			// advance ptr[d] past already assigned vertices
			while (ptr[d] < L.size() && !sol.available[L[ptr[d]]->index]) ++ptr[d];

			for (size_t k = ptr[d]; k < L.size() && examined < K; ++k) 
			{
				Ins::Vertex* cand = L[k];
				if (!sol.available[cand->index]) continue;

				double ndt; int brk;
				if (!try_append(d, cand, ndt, brk)) continue;

				// key: score / (Δt + capacity fraction)
				double delta_t = ndt - tour.deptime.back();
				double cap_frac = cand->weight / ins->t[d].W_max + cand->volume / ins->t[d].V_max;
				double denom = std::max(1e-12, delta_t) + std::max(1e-12, cap_frac);
				double key = cand->score / denom;

				if (key > best_key) 
				{
					best_key = key; best_v = cand; best_deptime = ndt; best_brk = brk; 
				}
				++examined;
			}

			if (best_v) 
			{
				// commit the seed
				tour.seq.push_back(best_v);
				tour.deptime.push_back(best_deptime);
				tour.max_shift.push_back(0);
				tour.action.push_back(best_brk);
				if (best_brk == 1) tour.breakindex = int(tour.seq.size()) - 1;
				tour.weight += best_v->weight;
				tour.volume += best_v->volume;
				tour.score += best_v->score;
				sol.score += best_v->score;
				sol.available[best_v->index] = false;
				seeded_any = true;
			}
		}
	}

	// 2) Phase B: greedily append remaining vertices
	std::vector<Ins::Vertex*> pool;
	pool.reserve(ins->maxvertices);
	for (int i = 1; i < ins->maxvertices - 1; ++i)
		if (sol.available[i]) pool.push_back(&ins->v[i]);
	std::sort(pool.begin(), pool.end(), [](Ins::Vertex* a, Ins::Vertex* b) { return a->score > b->score; });

	for (Ins::Vertex* cand : pool) {
		int best_d = -1; double best_key = -DBL_MAX;
		double best_deptime = 0.0; int best_brk = 0;

		for (int d = 0; d < ins->maxtours; ++d) 
		{
			double ndt; int brk;
			if (!try_append(d, cand, ndt, brk)) continue;

			double delta_t = ndt - sol.tours[d].deptime.back();
			double cap_frac = cand->weight / ins->t[d].W_max + cand->volume / ins->t[d].V_max;
			double denom = std::max(1e-12, delta_t) + std::max(1e-12, cap_frac);
			double key = cand->score / denom;

			if (key > best_key) { best_key = key; best_d = d; best_deptime = ndt; best_brk = brk; }
		}

		if (best_d >= 0) 
		{
			auto& tour = sol.tours[best_d];
			tour.seq.push_back(cand);
			tour.deptime.push_back(best_deptime);
			tour.max_shift.push_back(0);
			tour.action.push_back(best_brk);
			if (best_brk == 1) tour.breakindex = int(tour.seq.size()) - 1;
			tour.weight += cand->weight;
			tour.volume += cand->volume;
			tour.score += cand->score;
			sol.score += cand->score;
			sol.available[cand->index] = false;
		}
	}

	// 3) Close tours with end depot
	for (int t = 0; t < ins->maxtours; ++t) 
	{
		Sol::Tour& tour = sol.tours[t];
		double currenttime = ins->t[t].EDT + tour.deptime.back();
		double arrivaltime = ins->arrival_time(tour.seq.back()->con[ins->maxvertices - 1], currenttime);
		tour.deptime.push_back(arrivaltime - ins->t[t].EDT);
		tour.seq.push_back(&ins->v[ins->maxvertices - 1]);
		tour.max_shift.push_back(0);
		if (tour.breakindex == -1) {
			tour.action.push_back(1);
			tour.deptime.back() += ins->breakdur;
			tour.breakindex = int(tour.seq.size()) - 1;
		}
		else {
			tour.action.push_back(0);
		}
		tour.calc_maxshift();
	}
	cout << "start solution with score: " << sol.score << endl;
	//sol.check();
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

Res Tabu::solve(int bestknown)
{
	using clock = std::chrono::steady_clock;
	auto t0 = clock::now();
	s.reset();
	//s.read_from_file();
	//s.write_to_file();
	parallel_construct(s);
	gb = s;//set global best to initial solution
	ElitePool elites;
	elites.consider(s, ins->maxvertices);

	//criteria selector based on constraint pressure
	enum Crit { SCORE = 0, TIME = 1, VOLUME = 2, WEIGHT = 3, N_CRIT = 4 };
	struct CritStats { double ema_gain = 0.0; long used = 0; };
	std::array<CritStats, N_CRIT> crit_stats{};
	const double EMA_RHO = 0.1; //	learning rate for exponential moving average of gains
	const double THRESH = 0.80;   // start biasing after 80% utilization
	const double GAMMA = 5;    // bias strength
	const double EPS_GAIN = 1e-6;// to avoid zero gains

	uniform_int_distribution<> nbpicker(1,3);//random move selector
	TabuVector tabulist(ins->maxtours,ins->maxvertices,nb_tabu_it);//tabulist init

	int noimpr = 0;//number of iterations with no improvement
	int iter = 0;//number of iterations
	int nonb1 = 0;//non improving counter for move 1
	int nonb2 = 0;//non improving counter for move 2
	int nonb3 = 0;//non improving counter for move 3
	int no_feasible_moves_in_a_row = 0;
	int restart_count = 0;
	double restart_dist_sum = 0.0;


	while (noimpr < max_noimpr)
	{
		if (gb.score == ins->maxscore)
		{
			break;
		}
		Sol backup = s;;//backup current solution
		int prev_score = s.score;
		bool moved = false;
		//select criterion based on constraint pressure
		auto P = compute_pressures(s);
		auto ctx_mult = [&](int crit)->double {
			if (crit == TIME)   return 1.0 + std::max(0.0, P.time_p - THRESH) * GAMMA;
			if (crit == VOLUME) return 1.0 + std::max(0.0, P.volume_p - THRESH) * GAMMA;
			if (crit == WEIGHT) return 1.0 + std::max(0.0, P.weight_p - THRESH) * GAMMA;
			return 1.0; // SCORE neutral
			};

		auto perf = [&](int crit)->double { return crit_stats[crit].ema_gain + EPS_GAIN; };

		std::array<double, N_CRIT> w = 
		{
			perf(SCORE) * ctx_mult(SCORE),
			perf(TIME) * ctx_mult(TIME),
			perf(VOLUME) * ctx_mult(VOLUME),
			perf(WEIGHT) * ctx_mult(WEIGHT)
		};

		// 1) Sanitize: set negatives/NaN/inf to 0
		for (double& wi : w) {
			if (!(wi > 0.0) || !std::isfinite(wi)) wi = 0.0;  // catches NaN/inf/neg
		}

		// 2) Normalize by max so the largest weight is 1.0 (keeps magnitudes healthy)
		double maxw = *std::max_element(w.begin(), w.end());
		if (maxw > 0.0) {
			for (double& wi : w) wi /= maxw;
		}

		// 3) Mix with a uniform floor so nothing collapses to zero
		//    (p' = 0.85*p + 0.15*1.0)
		for (double& wi : w) wi = 0.85 * wi + 0.15;

		auto pick_weighted = [&](const std::array<double, N_CRIT>& ww)->int {
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
			return N_CRIT - 1; // fallback
			};

		int crit_id = pick_weighted(w);
		auto kind = static_cast<RatioKind>(crit_id);
		//select neighborhood structure at random
		int pick=nbpicker(engine);
		//int pick = 2;//to debug
		//build admissable neighborhoods using the selected neighborhoodstructure
		switch (pick)
		{
			case 1:
			{
				auto nb = one_one_replace_gen_nb(s,tabulist,gb.score,kind);
				if (nb.size() == 0)
				{
					//perturbe(s);
					//s.check();
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
				auto nb = two_one_replace_gen_nb(s,tabulist,gb.score,kind);
				if (nb.size() == 0)
				{
					//perturbe(s);
					//s.check();
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
				auto nb = one_two_replace_gen_nb(s,tabulist,gb.score,kind);
				if (nb.size() == 0)
				{
					//perturbe(s);
					//s.check();
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
		}//end switch
		/*
		if (!s.check())
		{
			cout << "error in replace: " << pick << endl;
		}
		*/
		if (no_feasible_moves_in_a_row > 5) 
		{
			int idx = elites.pick_idx(s, ins->maxvertices, engine);
			//int idx = elites.pick_farthest_idx(s, ins->maxvertices);
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
		
		// --- RVND: randomize LS order, restart when any op improves ---
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
					case LsOp::TwoOpt:   improved = two_opt_nb(s, 1);    break;   // or moves.two_opt_nb(s,1)
					case LsOp::Swap:     improved = swap_nb(s, 1);       break;
					case LsOp::Swap2:    improved = swap2_nb(s, 1);      break;
					case LsOp::Relocate: improved = relocate_nb(s, 1);   break;
				}
				if (improved) 
				{                     // RVND "restart-on-improvement"
					improved_any = true;
					break;
				}
			}
		}
		//update constraint pressure
		//cout << "it: " << iter << " score: " << s.score << " best: " << gb.score << " noimpr: " << noimpr << endl;
		if (moved) 
		{
			double gain = std::max(0, s.score - prev_score);
			crit_stats[crit_id].ema_gain =
				(1.0 - EMA_RHO) * crit_stats[crit_id].ema_gain + EMA_RHO * gain;
			crit_stats[crit_id].used++;
		}
		else 
		{
			// gentle global decay
			for (int c = 0; c < N_CRIT; ++c)
				crit_stats[c].ema_gain *= (1.0 - 0.05 * EMA_RHO);
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
		/*
		if ((iter % 2000) == 0 && restart_count > 0) 
		{
			std::cout<< s.score  << " [Restart] count=" << restart_count<< " avg_dist=" << (restart_dist_sum / restart_count) << "\n";
		}
		*/
		//cout << iter << endl;
	}//end while smaller than max_noimpr
	double cpuTime = std::chrono::duration<double>(clock::now() - t0).count();
	gb.check();
	//gb.write_to_cplex();
	std::cout<<"iter without replacement nb: " << nonb1<<" <> " << nonb2<<" <> " << nonb3 <<" total iterations: "<< iter << endl;
	//std::cout << gb << endl;
	return Res(gb, cpuTime, bestknown);
}


