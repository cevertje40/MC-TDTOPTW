#include "solver.h"

void Moves::calculate_maxshift(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int size = (int) sol.solution[d].size();
		double departuretime = 0;
		double arrivaltime = ins->t[d].LAT;
		if (sol.action[d].back() == 1)//break at enddepot
		{
			//traveltime contains breakdur
			sol.max_shift[d].back() = (ins->t[d].T_max - (sol.traveltime[d].back()-ins->breakdur));
		}
		else
		{
			sol.max_shift[d].back() = (ins->t[d].T_max - sol.traveltime[d].back());
		}
		Ins::Vertex* y;
		Ins::Vertex* z;
		for (int i = 0; i < size - 2; ++i)// depots don't count
		{
			//define the 2 elements
			y = sol.solution[d][size - (i + 2)];
			int breaki = sol.action[d][size - (i + 2)];
			z = sol.solution[d][size - (i + 1)];
			departuretime = ins->departure_time(y->con[z->index], arrivaltime);
			//TW check on departuretime 
			if (departuretime > y->UTW[d] + y->serv)
			{
				departuretime = y->UTW[d] + y->serv;
			}
			//break may limit the maximum allowable shift
			if (breaki == 1)
			{
				if (departuretime > ins->breakend + ins->breakdur)
				{
					cout<<"bij calc maxshift break verhindert een maxshift: " << departuretime <<"<=>"<< ins->breakend + ins->breakdur << endl;
					departuretime = ins->breakend + ins->breakdur;
				}
			}
			//store result
			sol.max_shift[d][size - (i + 2)] = departuretime - (sol.traveltime[d][size - (i + 2)] + ins->t[d].EDT);
			//reset variable for the calculation of next point
			arrivaltime = departuretime - (y->serv+breaki*ins->breakdur);
		}// end for i
	}//end for all trucks
}//end calc maxshift

void Moves::best_insert_nb(Sol& sol)
{
	//1.repeat until no improvement can be found
	bool improvement = true;
	int counter = 0;
	while (improvement)
	{
		improvement = false;
		double ratio = 0.0;
		int position = -1;
		Ins::Vertex* candidate=NULL;
		int bestpath = -1;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			//sort possible candidates
			int endj = (int)sol.solution[d].size();
			for (int j = 0; j < endj - 1; ++j)// for all  inlcuded vertices in the solution (non-depot)
			{
				Ins::Vertex* x = sol.solution[d][j];//point before insertion
				int nb_size = (int) x->nb[d].size();
				for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex
				{
					Ins::Vertex* y = x->nb[d][i];//point that might be inserted
					Ins::Vertex* z = sol.solution[d][j + 1];//point to shift
					int breakz = sol.action[d][j + 1];
					if ((sol.available[y->index]) && (y->nbi[d][z->index]))//y moet buur van z zijn want 
					{
						if ((sol.volume[d] + y->volume < ins->t[d].V_max) && (sol.weight[d] + y->weight < ins->t[d].W_max))
						{//increase in capacity moet nog gecontroleerd worden
							//gather departure time
							double currenttime = sol.traveltime[d][j] + ins->t[d].EDT;//service bij x zit hier al in
							//travel time from x to y
							double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
							if (arrivaltime < y->LTW[d])//break inserten kan niet dus break kan ltw niet dichter brengen
							{
								arrivaltime = y->LTW[d];
							}
							if (arrivaltime > y->UTW[d])
							{
								continue;//mag je al stoppen met rekenen voor dit punt op deze positie
							}
							arrivaltime += y->serv;//break inserten kan niet
							//travel time from y to z
							double waitz = 0;//wait at z
							arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
							if (arrivaltime + breakz * (ins->breakdur) < z->LTW[d])
							{
								waitz = z->LTW[d] - (arrivaltime + breakz * (ins->breakdur));
								arrivaltime = z->LTW[d] - breakz * (ins->breakdur);

							}
							arrivaltime += z->serv;
							double shift = (arrivaltime - ins->t[d].EDT) - sol.traveltime[d][j + 1];//increase in travel time
							if (shift <= sol.max_shift[d][j + 1])//check of het punt geinsert kan worden
							{
								improvement = true;
								double weightavail = ins->t[d].W_max - sol.weight[d];
								double volumeavail = ins->t[d].V_max - sol.volume[d];
								double ratiocheck;
								int function = 0;
								if (shift <= waitz)//als shift=0 of kleiner dan de weight is: speciaal regime toepassen
								{
									shift = 1;
									function = 1;
								}
								//ratiocheck = (double(y->score) / shift);
								//ratiocheck = (double(y->score*y->score) / shift);
								ratiocheck = (double(y->score * y->score) / shift) * (1 + function);
								//ratiocheck = (double(y->score*y->score) / shift)*(1 + function)*(weightavail / y->weight);
								//ratiocheck = (double(y->score*y->score) / shift)*(1 + function)*(volumeavail / y->volume);
								if (ratiocheck > ratio)//enkel op minimale increase checken
								{//update candidates
									ratio = ratiocheck;
									position = j;
									candidate = y;
									bestpath = d;
								}
							}// end if check feasible insertion
						}//if cap restrictions
					}//end if available
				}//end i
			}// end for j
		}//end for al paths
		if (improvement)
		{
			// execute best insertion if improvement found
			sol.available[candidate->index] = false;
			sol.solution[bestpath].insert(sol.solution[bestpath].begin() + position + 1, candidate);//insert point y after x
			sol.traveltime[bestpath].insert(sol.traveltime[bestpath].begin() + position + 1, 0);//insert temporary value
			sol.action[bestpath].insert(sol.action[bestpath].begin() + position + 1, 0);//insert regular visit action change later when necessary
			sol.scores[bestpath] += candidate->score;// update score of the new solution
			sol.score += candidate->score;// update score of the new solution
			sol.volume[bestpath] += candidate->volume;
			sol.weight[bestpath] += candidate->weight;
			if (position < sol.breakindex[bestpath])
			{
				sol.breakindex[bestpath] += 1;//update breakindex
			}
			//update travel time solution vector 
			sol.max_shift[bestpath].insert(sol.max_shift[bestpath].begin() + position + 1, 0);
			double currenttime = sol.traveltime[bestpath][position] + ins->t[bestpath].EDT;
			for (int u = position + 1; u < sol.solution[bestpath].size(); ++u)
			{
				//gather departure time and corresponding time slot
				Ins::Vertex* o = sol.solution[bestpath][u - 1];
				Ins::Vertex* p = sol.solution[bestpath][u];
				//travel time from van o to p
				double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
				if (arrivaltime + sol.action[bestpath][u] * (ins->breakdur) < p->LTW[bestpath])
				{
					arrivaltime = p->LTW[bestpath] - sol.action[bestpath][u] * (ins->breakdur);
				}
				arrivaltime += p->serv+sol.action[bestpath][u]*ins->breakdur;
				sol.max_shift[bestpath][u] = (sol.traveltime[bestpath][u] + sol.max_shift[bestpath][u]) - (arrivaltime - ins->t[bestpath].EDT);
				//old arrival time + old maxshift - new arrival time = new max_shift
				sol.traveltime[bestpath][u] = arrivaltime - ins->t[bestpath].EDT;
				currenttime = arrivaltime;
			}//end for
			//the value of max_shift j+1 is now wrong but will be soon be updated
			//complete update of max_shift for the included vertices before the insertion
			double departuretime = 0;
			double arrivaltime = (sol.traveltime[bestpath][position + 2] + ins->t[bestpath].EDT + sol.max_shift[bestpath][position + 2]) - sol.solution[bestpath][position + 2]->serv+sol.action[bestpath][position + 2]*ins->breakdur;//service time eraftrekken
			for (int vz = position + 1; vz > 0; --vz)
			{
				//define the 2 elements
				Ins::Vertex* f = sol.solution[bestpath][vz];
				int breakf = sol.action[bestpath][vz];
				Ins::Vertex* g = sol.solution[bestpath][vz + 1];
				//find proper departure time and time slot
				departuretime = ins->departure_time(f->con[g->index], arrivaltime);	//utw check
				if (departuretime > f->UTW[bestpath] + f->serv+breakf*ins->breakdur)
				{
					departuretime = f->UTW[bestpath] + f->serv+breakf*ins->breakdur;
				}
				if (breakf == 1)
				{
					if (departuretime > ins->breakend + ins->breakdur)
					{
						cout << "bij insert break verhindert een maxshift:  " << departuretime << "<=>" << ins->breakend + ins->breakdur << endl;
						departuretime = ins->breakend + ins->breakdur;
						//breakrepositon = true;//mogelijk kan de break vervroegd worden en zo moet de break maxshift niet afremmen
					}
				}
				sol.max_shift[bestpath][vz] = departuretime - (sol.traveltime[bestpath][vz] + ins->t[bestpath].EDT);
				departuretime -= f->serv;
				arrivaltime = departuretime;
			}// end for vz
			//try to pull break if the break is still positioned at the end depot
			int end = (int) sol.solution[bestpath].size();
			if ((sol.breakindex[bestpath] == end - 1) && (ins->t[bestpath].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "pull break" << endl;
				pull_break(sol, bestpath);
			}
		}//end improvement
	}//end while improvement
}

void Moves::pull_break(Sol& sol, int tour)
{
	int end = (int) sol.solution[tour].size();
	int earliestbreakindex = end - 1;
	for (int i = 0; i < end; ++i)
	{
		Ins::Vertex* p = sol.solution[tour][i];
		int breakcurrent = sol.action[tour][i];
		if ((ins->t[tour].EDT + sol.traveltime[tour][i]) - (p->serv+breakcurrent*ins->breakdur) >= ins->breakstart)
		{
			earliestbreakindex = i;
			break;
		}
	}
	for (int i = earliestbreakindex; i < sol.breakindex[tour]; ++i)//huidige breakindex niet evalueren die ken je al
	{
		double currenttime = ins->t[tour].EDT + sol.traveltime[tour][i - 1];
		bool feasible = true;
		for (int j = i - 1; j < end - 1; ++j)
		{
			Ins::Vertex* last = sol.solution[tour][j];
			Ins::Vertex* current = sol.solution[tour][j + 1];
			int breakcurrent = 0;
			if (j == i)
			{
				breakcurrent = 1;
			}
			double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
			if (arrivaltime + breakcurrent * (ins->breakdur) < current->LTW[tour])
			{
				arrivaltime = current->LTW[tour] - breakcurrent * (ins->breakdur);
			}
			//utw check want je kan later aankomen door de break vroeger te schedulen
			if (arrivaltime > current->UTW[tour])
			{
				feasible = false;
				break;
			}
			arrivaltime += current->serv+breakcurrent*ins->breakdur;
			currenttime = arrivaltime;
		}
		if ((feasible) && (currenttime <= ins->t[tour].EDT + sol.traveltime[tour].back()))
		{
			sol.action[tour][sol.breakindex[tour]] = 0;
			sol.action[tour][i] = 1;
			sol.breakindex[tour] = i;
			currenttime = ins->t[tour].EDT + sol.traveltime[tour][i - 1];
			for (int j = i - 1; j < end - 1; ++j)
			{
				Ins::Vertex* last = sol.solution[tour][j];
				Ins::Vertex* current = sol.solution[tour][j + 1];
				int breakcurrent = sol.action[tour][j + 1];
				double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
				if (arrivaltime + breakcurrent * (ins->breakdur) < current->LTW[tour])
				{
					arrivaltime = current->LTW[tour] - breakcurrent * (ins->breakdur);
				}
				arrivaltime += current->serv+breakcurrent*ins->breakdur;
				sol.traveltime[tour][j + 1] = arrivaltime - ins->t[tour].EDT;
				currenttime = arrivaltime;
			}
			sol.max_shift[tour].back() = (ins->t[tour].T_max - sol.traveltime[tour].back());
			double departuretime = 0;
			double arrivaltime = ins->t[tour].LAT;
			if (sol.action[tour].back() == 1)//break op enddepot
			{
				if (ins->t[tour].LAT > ins->breakend + ins->breakdur)
				{
					//cout<<"path: "<<d<< "pull break break op enddepot verhindert een maxshift: " << endl;
					sol.max_shift[tour].back() = (ins->breakend + ins->breakdur) - (sol.traveltime[tour].back() + ins->t[tour].EDT);
					arrivaltime = ins->breakend;//zoals hieronder service of enkel break in dit geval ervan aftrekken
				}
			}
			Ins::Vertex* y;
			Ins::Vertex* z;
			for (int i = 0; i < end - 2; ++i)// depots don't count
			{
				//define the 2 elements
				y = sol.solution[tour][end - (i + 2)];
				int breaki = sol.action[tour][end - (i + 2)];
				z = sol.solution[tour][end - (i + 1)];
				departuretime = ins->departure_time(y->con[z->index], arrivaltime);
				//TW check on departuretime 
				if (departuretime > y->UTW[tour] + y->serv+breaki*ins->breakdur)
				{
					departuretime = y->UTW[tour] + y->serv+breaki*ins->breakdur;
				}
				//break may limit the maximum allowable shift
				if (breaki == 1)
				{
					if (departuretime > ins->breakend + ins->breakdur)
					{
						//cout << "bij pull break break verhindert een maxshift: " << departuretime << "<=>" << breakend + breaktime << endl;
						departuretime = ins->breakend + ins->breakdur;
					}
				}
				//store result
				sol.max_shift[tour][end - (i + 2)] = departuretime - (sol.traveltime[tour][end - (i + 2)] + ins->t[tour].EDT);
				//reset variable for the calculation of next point
				arrivaltime = departuretime - (y->serv+breaki*ins->breakdur);
			}// end for i
			break;
		}//end for improvement
	}//end for all break combinations
}




Aco::Aco(Ins& ins, double alpha, double beta, double rho, int max_ants, int max_sol): Moves(ins), alpha(alpha), beta(beta), rho(rho), max_ants(max_ants)
{
	max_it = int(max_sol / max_ants);
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
				eta[i][j] = ins.v[j].score / (ins.v[i].con[j]->determin);
			}
			else
				eta[i][j] = 0.0;
		}
	}
	//creation of solutions
	gbs = Sol(ins);//ant class object
	ss.resize(max_ants);//array of ant class objects
	#pragma omp parallel for// parallel first touch to increase speed
	for (int i = 0; i < max_ants; ++i) 
	{
		ss[i] = Sol(ins);
	}// end for all ants
	iter_nr = 0;
	iter_score = -1;
}

void Aco::construct(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		bool breaktaken = false;
		while (sol.solution[d].back()->index != ins->maxvertices - 1) //until one solution is full=>sequential procedure
		{
			vector<double> prob_v;
			vector<double> temp_traveltime;
			vector<int> temp_action;
			Ins::Vertex* last = sol.solution[d].back();
			prob_v.resize(ins->maxvertices);
			temp_traveltime.resize(ins->maxvertices);
			temp_action.resize(ins->maxvertices);
			for (int i = 0; i < ins->maxvertices; ++i)
			{//check 3 constraints for potential vertex to include: availability,weight,volume
				if (i != ins->maxvertices - 1)
				{
					if (sol.available[i])
					{
						if ((sol.volume[d] + ins->v[i].volume < ins->t[d].V_max) && (sol.weight[d] + ins->v[i].weight < ins->t[d].W_max))
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
					double currenttime = sol.traveltime[d].back() + ins->t[d].EDT;
					Ins::Vertex* neighbor = &ins->v[i];
					double arrivaltime = ins->arrival_time(last->con[neighbor->index], currenttime);
					int brk = -1;
					//als je na breakstart aankomt moet je breaken
					if ((breaktaken == false) && (arrivaltime >= ins->breakstart))
					{
						arrivaltime += ins->breakdur;
						brk = 1;
					}
					else
					{
						brk = 0;
					}
					if (arrivaltime < neighbor->LTW[d])
					{
						prob_v[i] = 1 - (double(neighbor->LTW[d] - arrivaltime) / ins->t[d].T_max);//wachttijd meegeven
						arrivaltime = neighbor->LTW[d];//wachten als je te vroeg bent	
					}
					if (arrivaltime > neighbor->UTW[d])
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
					if (enddepottime > ins->t[d].LAT)//enddepot heeft geen service time
					{
						prob_v[i] = 0;//discard vertex if infeasible
					}
					else
					{
						if (neighbor->index != ins->maxvertices - 1)
						{
							enddepot = false;//if another point than the enddepot is found
						}
						temp_traveltime[i] = arrivaltime - ins->t[d].EDT;
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
			sol.solution[d].push_back(&ins->v[sel - 1]);
			sol.traveltime[d].push_back(temp_traveltime[sel - 1]);
			sol.available[sel - 1] = false;
			sol.scores[d] += ins->v[sel - 1].score;
			sol.score += ins->v[sel - 1].score;
			sol.max_shift[d].push_back(0);//dummy die dan in calc max shift upgedate wordt
			sol.volume[d] += ins->v[sel - 1].volume;
			sol.weight[d] += ins->v[sel - 1].weight;
			sol.action[d].push_back(temp_action[sel - 1]);
			if (temp_action[sel - 1] == 1)
			{
				breaktaken = true;
				sol.breakindex[d] = int(sol.solution[d].size()) - 1;
			}
			prob_v[sel - 1] = 0;
			//++counter;
		}// end while

		sol.available[ins->maxvertices - 1] = true;// end depot moet terug available zijn voor de volgende tour
		if (breaktaken == false)
		{
			sol.action[d].back() = 1;
			sol.traveltime[d].back() += ins->breakdur;
			sol.breakindex[d] = int(sol.solution[d].size()) - 1;
		}
	}//end for all d
	//end depot is now unavailable for inserts and replacements
	sol.available[ins->maxvertices - 1] = false;
}

void Aco::solve()
{
	construct(ss[0]);
	swap_nb(ss[0]);
	ss[0].check();
	calculate_maxshift(ss[0]);
	best_insert_nb(ss[0]);
	ss[0].check();
	cout << "debug" << endl;

}

void Moves::swap_nb(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		bool improvement = true;
		while (improvement)
		{
			improvement = false;
			int end = (int)sol.solution[d].size();
			for (int i = 1; i < end - 2; ++i)// every vertex except depots and second last real vertex
			{
				for (int j = i + 1; j < end - 1; ++j)// every vertex except depots and second last real vertex
				{
					Ins::Vertex* k = sol.solution[d][i - 1];//pred swap partner 1
					Ins::Vertex* y = sol.solution[d][i];//swap partner 1
					Ins::Vertex* l = sol.solution[d][i + 1];//suc partner 1
					Ins::Vertex* m = sol.solution[d][j - 1];//pred swap partner 2
					Ins::Vertex* z = sol.solution[d][j];//swap partner 2
					Ins::Vertex* n = sol.solution[d][j + 1];//suc partner 2
					if (z == l)
					{
						l = sol.solution[d][j - 1];//suc partner 1
						m = sol.solution[d][i + 1];//pred swap partner 2
					}
					if ((k->nbi[d][z->index]) && (m->nbi[d][y->index]) && (z->nbi[d][l->index]) && (y->nbi[d][n->index]))
					{//geen cap restricties noding hier blijft gelijk
						double delta_tt = 0;
						double departuretime = ins->t[d].EDT + sol.traveltime[d][i - 1];//bij punt voor y
						double currenttime = departuretime;
						double newtraveltime;
						double oldtraveltime;
						//tijdelijk swappen
						sol.solution[d][i] = z;
						sol.solution[d][j] = y;
						bool reqbreak = false;
						if (i <= sol.breakindex[d])
						{//break zit na i
							reqbreak = true;
						}
						for (int l = i - 1; l < j + 1; ++l)//local evaluation
						{
							Ins::Vertex* first = sol.solution[d][l];
							Ins::Vertex* second = sol.solution[d][l + 1];
							double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
							//do not account for break repositioning
							if (arrivaltime + sol.action[d][l + 1] * (ins->breakdur) < second->LTW[d])
							{
								arrivaltime = second->LTW[d] - sol.action[d][l + 1] * (ins->breakdur);
							}
							if (arrivaltime > second->UTW[d])//als hij tussen  z en y over een utw gaat moet je hem stoppen
							{
								goto stop;//swap is infeasible, stop for this combination
							}
							arrivaltime += second->serv+sol.action[d][l + 1]*ins->breakdur;
							currenttime = arrivaltime;
						}// end l
						newtraveltime = currenttime - departuretime;
						oldtraveltime = ((ins->t[d].EDT + sol.traveltime[d][j]) - departuretime);
						delta_tt = newtraveltime - oldtraveltime;
						if (delta_tt < 0)
						{
							currenttime = departuretime;
							for (int l = i - 1; l < end - 1; ++l)//global evaluation
							{
								Ins::Vertex* first = sol.solution[d][l];
								Ins::Vertex* second = sol.solution[d][l + 1];
								double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
								//account for break
								if ((reqbreak) && (arrivaltime >= ins->breakstart))
								{
									arrivaltime += ins->breakdur;
									reqbreak = false;
								}
								if (arrivaltime < second->LTW[d])
								{
									arrivaltime = second->LTW[d];
								}
								if (arrivaltime > second->UTW[d])//als hij tussen  z en y over een utw gaat moet je hem stoppen
								{
									goto stop;//swap is infeasible, stop for this combination
								}
								arrivaltime += second->serv;
								currenttime = arrivaltime;
							}// end l
							if (reqbreak == true)//als je aankomt bij het einddepot en nog steeds geen break genomen hebt
							{
								currenttime += ins->breakdur;//breaktime bijtellen bij aankomst tijd bij einddepot
							}
							newtraveltime = currenttime - departuretime;
							oldtraveltime = ((ins->t[d].EDT + sol.traveltime[d][end - 1]) - departuretime);
							delta_tt = newtraveltime - oldtraveltime;
							if (delta_tt < 0)
							{
								improvement = true;
								//double remember=sol.traveltime[d].back();
								//cout<<"feasible swap possible with improvement of"<<delta_tt<<endl;
								//travel time herevalueren
								currenttime = departuretime;
								//reset break variable
								bool reqbreak = false;
								if (i <= sol.breakindex[d])
								{
									reqbreak = true;
								}
								for (int u = i; u < end; ++u)
								{
									Ins::Vertex* first = sol.solution[d][u - 1];
									Ins::Vertex* second = sol.solution[d][u];
									double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
									//add possible waiting time
									if ((reqbreak) && (arrivaltime >= ins->breakstart))
									{
										arrivaltime += ins->breakdur;
										sol.action[d][u] = 1;
										sol.breakindex[d] = u;
										reqbreak = false;
									}
									else
									{
										sol.action[d][u] = 0;
									}
									if (arrivaltime < second->LTW[d])
									{
										arrivaltime = second->LTW[d];
									}
									arrivaltime += second->serv;
									sol.traveltime[d][u] = arrivaltime - ins->t[d].EDT;
									currenttime = arrivaltime;
								}//end for u
								//if by shortening the route no break is necessary anymore, set break at enddepot
								if (reqbreak)
								{
									sol.action[d].back() = 1;
									sol.traveltime[d].back() += ins->breakdur;
									sol.breakindex[d] = end - 1;
								}
								//cout<<"actual improvement"<<sol.traveltime[d].back()-remember<<endl;
								//check_solution(sol);
								//cout << "debug" << endl;
								break;
								//i is geswapped ofwel j lus terug doorlopen ofwel naar volgende i gaan
							}//end executed swap
							else
							{
								//terugswappen als het niet gaat
								sol.solution[d][i] = y;
								sol.solution[d][j] = z;
							}
						}
						else
						{
						stop:
							//terugswappen als het niet gaat
							sol.solution[d][i] = y;
							sol.solution[d][j] = z;
						}
					}
				}//end j
			}// end i
		}//end while improvement
	}//end for all d

}
