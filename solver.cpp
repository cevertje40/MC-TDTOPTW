#include "solver.h"

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
				eta[i][j] = max(0.001,ins.v[j].score) / consumption;
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
	output << "score: " << res.score << " after: " << res.time << " gap: " << res.gap << "\n";
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
				move_nb(s[ant], 1);
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
	cout << gb << endl;
	//cout << "best score: "<<gb.score<<"after: "<<cpuTime << endl;
	return Res(gb.score, cpuTime, bestknown);
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
					double consumption = ((last->con[i]->determin / ins->t[tour.index].T_max) + (ins->v[i].weight / ins->t[tour.index].W_max) + (ins->v[i].volume / ins->t[tour.index].V_max)) / 3;
					prob_v[i] = max(0.001,ins->v[i].score) / consumption;
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
		move_nb(s,1);
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
	cout << gb << endl;
	//cout << "best score: "<<gb.score<<"after: "<<cpuTime << endl;
	return Res(gb.score, cpuTime, bestknown);
}


Tabu::Tabu(Ins& ins, int max_noimpr, int nb_tabu_it): Moves(ins),max_noimpr(max_noimpr), nb_tabu_it(nb_tabu_it)
{
	gb = Sol(ins);//best sol
	s = Sol(ins);//iter sol
}

void Tabu::perturbe(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		//find most consuming vertex
		int targetvertex = -1;
		int end = (int) sol.tours[d].seq.size();
		double maxratio= -DBL_MAX;
		for (int i = 1; i < end-1; ++i)
		{
			Ins::Vertex* a = sol.tours[d].seq[i - 1];
			Ins::Vertex* b = sol.tours[d].seq[i];
			double usedtraveltime = (sol.tours[d].deptime[i]-(b->serv+sol.tours[d].action[i]*ins->breakdur)) - sol.tours[d].deptime[i - 1];//includes waiting time
			double mintraveltime = a->con[b->index]->determin;
			double ratio = usedtraveltime / mintraveltime;
			if (ratio > maxratio)
			{
				targetvertex = i;
				maxratio = ratio;
			}
		}
		//delete most consuming verter
		if (targetvertex != -1)
		{
			sol.remove_vertex(sol.tours[d], targetvertex);
		}
	}
}

void Tabu::parallel_construct(Sol& sol)
{
	//only use candidates that you can reach from the stard depot and still return to the end depot
	vector<Ins::Vertex*>candidates;
	for (int i = 1; i < ins->maxvertices-1; ++i)
	{
		if (ins->v[0].con[i]->determin+ins->v[i].serv+ins->v[i].con[ins->maxvertices-1]->determin <= ins->t[0].T_max)
		{
			candidates.push_back(&ins->v[i]);
		}
	}
	for (int i = 0; i < candidates.size() - 1; ++i)
	{
		//calculate traveltime for candidate under consideration
		vector<double>traveltime(ins->maxtours, 0.0);//traveltime per route
		vector<bool>feasible(ins->maxtours, false);//feasibility per route
		int besttour = -1;
		double besttraveltime = DBL_MAX;
		int bestbrk = 0;
		bool succes = false;
		for (int t = 0; t < ins->maxtours; ++t)
		{
			Sol::Tour& tour = sol.tours[t];
			double currenttime = tour.deptime.back() + ins->t[t].EDT;
			Ins::Vertex* last = tour.seq.back();
			Ins::Vertex* candidate = candidates[i];
			if ((tour.weight + candidate->weight < ins->t[t].W_max) && (tour.volume + candidate->volume < ins->t[t].V_max))
			{
				double arrivaltime = ins->arrival_time(last->con[candidate->index], currenttime);
				int brk = -1;
				//check if break is needed
				if ((tour.breakindex == -1) && ((arrivaltime >= ins->breakstart) || (candidate->LTW[t] - ins->breakdur >= ins->breakstart)))
				{
					arrivaltime += ins->breakdur;
					brk = 1;
				}
				else
				{
					brk = 0;
				}
				//tw checks
				if (arrivaltime < candidate->LTW[t])
				{
					arrivaltime = candidate->LTW[t];//wachten als je te vroeg bent	
				}
				if (arrivaltime > candidate->UTW[t])
				{
					continue;//stop the calculation
					feasible[t] = false;
				}
				//add service time
				arrivaltime += candidate->serv;
				//return to end depot check
				double enddepottime = ins->arrival_time(candidate->con[ins->maxvertices - 1], arrivaltime);
				if ((tour.breakindex == -1) && (brk == 0))//no break taken yet and you are not going to break at candidate
				{
					enddepottime += ins->breakdur;//take break at end depot
				}
				if (enddepottime > ins->t[t].LAT)//enddepot heeft geen service time
				{
					feasible[t] = false;
				}
				else
				{
					succes = true;
					feasible[t] = true;
					traveltime[t] = ((arrivaltime - brk * ins->breakdur) - currenttime);//exclude break as this would be unfair when no break is necessary for some candidates
					if (traveltime[t] < besttraveltime)//candidates are already sorted from high score to low score
					{
						besttour = t;
						besttraveltime = arrivaltime - ins->t[t].EDT;//include the possible break time
						bestbrk = brk;
					}
				}//end else
			}//end cap constraints
		}//end for all tours
		//insert candidate with highest score in a route with the lowest traveltime increase
		if (succes)
		{
			Sol::Tour& tour = sol.tours[besttour];
			tour.seq.push_back(candidates[i]);
			tour.deptime.push_back(besttraveltime);
			sol.available[candidates[i]->index] = false;
			tour.score += candidates[i]->score;
			sol.score += candidates[i]->score;
			tour.max_shift.push_back(0);//dummy die dan in calc max shift upgedate wordt
			tour.volume += candidates[i]->volume;
			tour.weight += candidates[i]->weight;
			tour.action.push_back(bestbrk);
			if (bestbrk == 1)
			{
				tour.breakindex = int(tour.seq.size()) - 1;
			}
		}
	}
	//add end depot to all routes
	for (int t = 0; t < ins->maxtours; ++t)
	{
		Sol::Tour& tour = sol.tours[t];
		double currenttime = ins->t[t].EDT + tour.deptime.back();
		double arrivaltime = ins->arrival_time(tour.seq.back()->con[ins->maxvertices - 1], currenttime);
		tour.deptime.push_back(arrivaltime - ins->t[t].EDT);
		tour.seq.push_back(&ins->v[ins->maxvertices - 1]);
		tour.max_shift.push_back(0);//dummy die dan in calc max shift upgedate wordt
		if (tour.breakindex == -1)
		{
			tour.action.push_back(1);
			tour.deptime.back() += ins->breakdur;
			tour.breakindex = int(tour.seq.size()) - 1;
		}
		else
		{
			tour.action.push_back(0);
		}
		tour.calc_maxshift();
	}
	//sol.check();
}//end parallel construct

template<typename MoveType>
void executeMove(boost::heap::priority_queue<MoveType>& nb, Sol& s, TabuVector& tabulist) 
{
	if (nb.size() >= 1)
	{
		auto exec_nb = nb.top();
		exec_nb.execute(s);
		for (int i = 0; i < exec_nb.move.first.size(); ++i)
		{
			tabulist.addTabu(exec_nb.move.first[i]->index,exec_nb.tour);
		}
		for (int i = 0; i < exec_nb.move.second.size(); ++i)
		{
			tabulist.addTabu(exec_nb.move.second[i]->index,exec_nb.tour);
		}
		tabulist.nextIteration();
	}
}


Res Tabu::solve(int bestknown)
{
	clock_t start, end;
	start = clock();
	int noimpr = 0;
	s.reset();
	//s.read_from_file();
	//s.write_to_file();
	parallel_construct(s);
	int debug_iter = 0;
	uniform_int_distribution<> nbpicker(1,3);
	TabuVector tabulist(ins->maxvertices,nb_tabu_it);
	ratiofunctions.push_back(&Moves::score);
	ratiofunctions.push_back(&Moves::score_tt);
	ratiofunctions.push_back(&Moves::score_v);
	double alpha = 0.9;
	double beta = 0.05;
	double gamma = 0.05;
	int nonb1 = 0;
	int nonb2 = 0;
	int nonb3 = 0;
	while (noimpr < max_noimpr)
	{
		//select neighborhood structure at random
		int pick=nbpicker(engine);
		//int pick = 2;
		//double alpha= rand() / (double)RAND_MAX;
		//double beta = rand() / (double)RAND_MAX;
		//double gamma= rand() / (double)RAND_MAX;
		vector<ScoreFunctionPointer> out;
		sample(ratiofunctions.begin(),ratiofunctions.end(),back_inserter(out),1,engine);
		//build admissable neighborhoods using the selected neighborhoodstructure
		switch (pick)
		{
			case 1:
			{
				auto nb = one_one_replace_gen_nb(s,tabulist,gb.score);
				//auto nb = one_one_replace_gen_nb(s, 50, alpha, beta, gamma,out[0]);
				//tabulist.printTabuList();
				if (nb.size() == 0)
				{
					perturbe(s);
					++nonb1;
				}
				executeMove(nb, s, tabulist);
				break;
			}
			case 2:
			{
				auto nb = two_one_replace_gen_nb(s,tabulist,gb.score);
				//auto nb = two_one_replace_gen_nb(s,50,alpha, beta, gamma,out[0]);
				if (nb.size() == 0)
				{
					perturbe(s);
					++nonb2;
				}
				executeMove(nb, s, tabulist);
				break;
			}
			case 3:
			{
				auto nb = one_two_replace_gen_nb(s,tabulist,gb.score);
				//auto nb = one_two_replace_gen_nb(s, 50, alpha, beta, gamma,out[0]);
				if (nb.size() == 0)
				{
					perturbe(s);
					++nonb3;
				}
				executeMove(nb, s, tabulist);
				break;
			}
		}//end switch
		if (!s.check())
		{
			cout << "error in replace" << endl;
		}
		//VND
		two_opt_nb(s, 1);
		if (!s.check())
		{
			cout << "error in two-opt" << endl;
		}
		swap_nb(s, 1);
		if (!s.check())
		{
			cout << "error in swap_nb" << endl;
		}
		swap2_nb(s, 1);
		if (!s.check())
		{
			cout << "error in swap2_nb" << endl;
		}
		move_nb(s, 1);
		if (!s.check())
		{
			cout << "error in move_nb" << endl;
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
		++debug_iter;
		//cout << debug_iter << endl;
	}//end while smaller than max_noimpr
	end = clock();
	double cpuTime;
	cpuTime = difftime(end, start) / CLOCKS_PER_SEC;
	gb.check();
	gb.write_to_cplex();
	cout<<"it with no nb: " << nonb1<<" <> " << nonb2<<" <> " << nonb3 << endl;
	cout << gb << endl;
	return Res(gb.score, cpuTime, bestknown);
}


