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
				double consumption = ((ins.v[i].con[j]->determin / ins.t[0].T_max) + (ins.v[j].weight / ins.t[0].W_max) + (ins.v[j].volume / ins.t[0].V_max)) / 3;
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
		sol.calc_maxshift(tour);
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
				two_opt_nb(s[ant], 1);
				swap_nb(s[ant],1);
				//swap2_nb(s[ant], 1);
				//move_nb(s[ant], 1);
				insert_nb(s[ant], 1);
				replace_nb(s[ant], 1);
				two_opt_nb(s[ant], 1);
				swap_nb(s[ant], 1);
				//swap2_nb(s[ant], 1);
				//move_nb(s[ant], 1);
				insert_nb(s[ant], 1);
				replace_nb(s[ant], 1);
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

void Ils::construct(Sol& sol)
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
		sol.calc_maxshift(tour);
	}//end for all d
	//make end depot unaivailable for other moves
	sol.available[ins->maxvertices - 1] = false;
	//sol.check();
}

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
				for (int vv = 0; vv < tour.seq.size(); ++vv)
				{
					tour.action[vv] = 0;
				}
				sol.update_traveltime_break(tour.index, post - 1, int(tour.seq.size()));

			}
			else
			{
				sol.update_traveltime(tour.index, post - 1, int(tour.seq.size()));
			}
			sol.calc_maxshift(tour);
		}//end if tour is not empty
	}//end for all tours
}//end shake

Ils::Ils(Ins& ins, int max_sol, int threshold1, int threshold2, int threshold3): Moves(ins), threshold1(threshold1), threshold2(threshold2), threshold3(threshold3)
{
	max_it = max_sol;
	gb = Sol(ins);//best sol
	s = Sol(ins);//iter sol
	iter_nr = 0;
	iter_score = -1;
	
}

Res Ils::solve(int bestknown)
{
	clock_t start, end;
	start = clock();
	s.reset();
	construct(s);
	int noimpr = 0;
	int post = 1;//post -->position to start the removal
	int cons = 1;//cons -->amount of vertices to be removed
	for (int iter = 0; iter < max_it; ++iter)//iteration loop
	{
		if ((noimpr > threshold2) && ((noimpr + 1) % threshold3 == 0))
		{
			next_permutation(s.tourindex.begin(), s.tourindex.end());//permutates the path index that is stored in the solution
		}
		else
		{
			shake(s, post, cons);
		}

		two_opt_nb(s, 1);
		swap_nb(s, 1);
		swap2_nb(s, 1);
		move_nb(s, 1);
		insert_nb(s, 1);
		replace_nb(s, 1);
		/*
		int cont = 1;
		while (cont >= 1)
		{
			cont = 0;
			cont += two_opt_nb(s, 1);
			cont += swap_nb(s, 1);
			cont += swap2_nb(s, 1);
			cont += move_nb(s, 1);
			cont += insert_nb(s, 1);
			cont += replace_nb(s, 1);
			//cout << "debug here" << endl;
		}
		*/
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
