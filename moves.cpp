#include "moves.h"

void Moves::calculate_maxshift(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int size = (int)sol.solution[d].size();
		double departuretime = 0;
		sol.max_shift[d].back() = (ins->t[d].T_max - sol.traveltime[d].back());
		double arrivaltime = ins->t[d].LAT - sol.action[d].back() * ins->breakdur;//if you break at the end depot subtract breakduration
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
					cout << "bij calc maxshift break verhindert een maxshift: " << departuretime << "<=>" << ins->breakend + ins->breakdur << endl;
					departuretime = ins->breakend + ins->breakdur;
				}
			}
			//store result
			sol.max_shift[d][size - (i + 2)] = departuretime - (sol.traveltime[d][size - (i + 2)] + ins->t[d].EDT);
			//reset variable for the calculation of next point
			arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
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
		Ins::Vertex* candidate = NULL;
		int bestpath = -1;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			//sort possible candidates
			int endj = (int)sol.solution[d].size();
			for (int j = 0; j < endj - 1; ++j)// for all  inlcuded vertices in the solution (non-depot)
			{
				Ins::Vertex* x = sol.solution[d][j];//point before insertion
				int nb_size = (int)x->nb[d].size();
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
								waitz = z->LTW[d] - (arrivaltime + breakz * ins->breakdur);
								arrivaltime = z->LTW[d] - breakz * (ins->breakdur);

							}
							arrivaltime += z->serv + breakz * ins->breakdur;
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
				sol.breakindex[bestpath] += 1;//due to insertion of 1 vertex the index needs to be incremented with 1
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
				if (arrivaltime + sol.action[bestpath][u] * ins->breakdur < p->LTW[bestpath])
				{
					arrivaltime = p->LTW[bestpath] - sol.action[bestpath][u] * ins->breakdur;
				}
				arrivaltime += p->serv + sol.action[bestpath][u] * ins->breakdur;
				sol.max_shift[bestpath][u] = (sol.traveltime[bestpath][u] + sol.max_shift[bestpath][u]) - (arrivaltime - ins->t[bestpath].EDT);
				//old arrival time + old maxshift - new arrival time = new max_shift
				sol.traveltime[bestpath][u] = arrivaltime - ins->t[bestpath].EDT;
				currenttime = arrivaltime;
			}//end for
			//the value of max_shift j+1 is now wrong but will be soon be updated
			//complete update of max_shift for the included vertices before the insertion
			double departuretime = 0;
			double arrivaltime = (sol.traveltime[bestpath][position + 2] + ins->t[bestpath].EDT + sol.max_shift[bestpath][position + 2]) - (sol.solution[bestpath][position + 2]->serv + sol.action[bestpath][position + 2] * ins->breakdur);//service time eraftrekken
			for (int vz = position + 1; vz > 0; --vz)
			{
				//define the 2 elements
				Ins::Vertex* f = sol.solution[bestpath][vz];
				int breakf = sol.action[bestpath][vz];
				Ins::Vertex* g = sol.solution[bestpath][vz + 1];
				//find proper departure time and time slot
				departuretime = ins->departure_time(f->con[g->index], arrivaltime);	//utw check
				if (departuretime > f->UTW[bestpath] + f->serv + breakf * ins->breakdur)
				{
					departuretime = f->UTW[bestpath] + f->serv + breakf * ins->breakdur;
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
			int end = (int)sol.solution[bestpath].size();
			if ((sol.breakindex[bestpath] == end - 1) && (ins->t[bestpath].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "pull break" << endl;
				pull_break(sol, bestpath);
			}
			sol.check();
			cout << "hier" << endl;
		}//end improvement
	}//end while improvement
}

void Moves::pull_break(Sol& sol, int tour)
{
	int end = (int)sol.solution[tour].size();
	int earliestbreakindex = end - 1;
	for (int i = 0; i < end; ++i)
	{
		Ins::Vertex* p = sol.solution[tour][i];
		int breakcurrent = sol.action[tour][i];
		if ((ins->t[tour].EDT + sol.traveltime[tour][i]) - (p->serv + breakcurrent * ins->breakdur) >= ins->breakstart)
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
				arrivaltime = current->LTW[tour] - (breakcurrent * ins->breakdur);
			}
			//utw check want je kan later aankomen door de break vroeger te schedulen
			if (arrivaltime > current->UTW[tour])
			{
				feasible = false;
				break;
			}
			arrivaltime += current->serv + breakcurrent * ins->breakdur;
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
					arrivaltime = current->LTW[tour] - (breakcurrent * ins->breakdur);
				}
				arrivaltime += current->serv + breakcurrent * ins->breakdur;
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
				if (departuretime > y->UTW[tour] + y->serv + breaki * ins->breakdur)
				{
					departuretime = y->UTW[tour] + y->serv + breaki * ins->breakdur;
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
				arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
			}// end for i
			break;
		}//end for improvement
	}//end for all break combinations
}