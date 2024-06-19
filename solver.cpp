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
				eta[i][j] = max(0.001,ins.v[j].score) / (ins.v[i].con[j]->determin);
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
						if ((tour.volume + ins->v[i].volume < ins->t[t].V_max) && (tour.weight + ins->v[i].weight < ins->t[t].W_max))
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
			double random = dist(mt);
			while (random >= total)
			{
				total += prob_v[sel];
				++sel;
			}
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
				insert_nb(s[ant], 1);
				replace_nb(s[ant], 1);
				swap2_nb(s[ant], 1);
				move_nb(s[ant], 1);
				insert_nb(s[ant], 1);
			}
		}
		update_global_best();//best solution is stored
		pheromone_update();// iteration best solution its arcs are augmented
	}
	end = clock();
	double cpuTime;
	cpuTime = difftime(end, start) / CLOCKS_PER_SEC;
	gb.check();
	cout << gb << endl;
	cout << "best score: "<<gb.score<<"after: "<<cpuTime << endl;
	return Res(gb.score,cpuTime,bestknown);
}

Ils::Ils(Ins& ins, int max_it,int max_ni): Moves(ins), max_it(max_it), max_ni(max_ni)
{

}
