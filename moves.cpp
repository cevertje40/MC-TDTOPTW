#include "moves.h"

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
			sol.insertvertex(bestpath, candidate, position);
			//try to pull break if the break is still positioned at the end depot
			int end = (int)sol.solution[bestpath].size();
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

void Moves::reschedule_breaks(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int end = (int) sol.solution[d].size();
		int earliestbreakindex = end - 1;
		int bestindex = -1;
		double bestdecrease = 0;
		bool improvement = false;
		for (int i = 0; i < end; ++i)
		{
			Ins::Vertex* p = sol.solution[d][i];
			int breakcurrent = sol.action[d][i];
			if ((ins->t[d].EDT + sol.traveltime[d][i]) - (p->serv+breakcurrent*ins->breakdur) >= ins->breakstart)
			{
				earliestbreakindex = i;
				break;
			}
		}
		for (int i = earliestbreakindex; i < sol.breakindex[d]; ++i)//huidige breakindex niet evalueren die ken je al
		{
			double currenttime = ins->t[d].EDT + sol.traveltime[d][i - 1];
			bool feasible = true;
			for (int j = i - 1; j < end - 1; ++j)
			{
				Ins::Vertex* last = sol.solution[d][j];
				Ins::Vertex* current = sol.solution[d][j + 1];
				int breakcurrent = 0;
				if ((j + 1) == i)
				{
					breakcurrent = 1;
				}
				double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
				if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[d])
				{
					arrivaltime = current->LTW[d] - (breakcurrent * ins->breakdur);
				}
				//utw check is nodig want je komt later aan
				if (arrivaltime > current->UTW[d])
				{
					feasible = false;
					break;
				}
				arrivaltime += current->serv+breakcurrent*ins->breakdur;
				currenttime = arrivaltime;
			}
			if ((feasible) && (currenttime < ins->t[d].EDT + sol.traveltime[d].back()))
			{
				double decrease = (ins->t[d].EDT + sol.traveltime[d].back()) - currenttime;
				if (decrease > bestdecrease)
				{
					bestindex = i;
					bestdecrease = decrease;
					improvement = true;
				}
			}//end for improvement
		}//end for all break combinations
		if (improvement)//execute best breakposition
		{
			sol.action[d][sol.breakindex[d]] = 0;
			sol.action[d][bestindex] = 1;
			sol.breakindex[d] = bestindex;
			double currenttime = ins->t[d].EDT + sol.traveltime[d][bestindex - 1];
			for (int j = bestindex - 1; j < end - 1; ++j)
			{
				Ins::Vertex* last = sol.solution[d][j];
				Ins::Vertex* current = sol.solution[d][j + 1];
				int breakcurrent = sol.action[d][j + 1];
				double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
				if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[d])
				{
					arrivaltime = current->LTW[d] - (breakcurrent * ins->breakdur);
				}
				arrivaltime += current->serv+breakcurrent*ins->breakdur;
				sol.traveltime[d][j + 1] = arrivaltime - ins->t[d].EDT;
				currenttime = arrivaltime;
			}
			sol.max_shift[d].back() = (ins->t[d].T_max - sol.traveltime[d].back());
			double departuretime = 0;
			double arrivaltime = ins->t[d].LAT;
			if (sol.action[d].back() == 1)//break op enddepot
			{
				if (ins->t[d].LAT > ins->breakend + ins->breakdur)
				{
					//cout<<"path: "<<d<< "pull break break op enddepot verhindert een maxshift: " << endl;
					sol.max_shift[d].back() = (ins->breakend + ins->breakdur) - (sol.traveltime[d].back() + ins->t[d].EDT);
					arrivaltime = ins->breakend;//zoals hieronder service of enkel break in dit geval ervan aftrekken
				}
			}
			Ins::Vertex* y;
			Ins::Vertex* z;
			for (int i = 0; i < end - 2; ++i)// depots don't count
			{
				//define the 2 elements
				y = sol.solution[d][end - (i + 2)];
				int breaki = sol.action[d][end - (i + 2)];
				z = sol.solution[d][end - (i + 1)];
				departuretime = ins->departure_time(y->con[z->index], arrivaltime);
				//TW check on departuretime 
				if (departuretime > y->UTW[d] + y->serv+breaki*ins->breakdur)
				{
					departuretime = y->UTW[d] + y->serv+breaki*ins->breakdur;
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
				sol.max_shift[d][end - (i + 2)] = departuretime - (sol.traveltime[d][end - (i + 2)] + ins->t[d].EDT);
				//reset variable for the calculation of next point
				arrivaltime = departuretime - y->serv+breaki*ins->breakdur;
			}// end for i
			sol.check();
			cout << "hier" << endl;
		}//end if improvement
	}//end for all tours
}//end reschedule_breaks

void Moves::best_replace_nb(Sol& sol)
{

}

void Moves::replace_nb(Sol& sol)
{
	//1.repeat until no improvement can be found
	bool improvement = true;
	while (improvement)
	{
		improvement = false;
		double ratio = 0.0;
		int position = -1;
		Ins::Vertex* candidate = NULL;
		double ttxybest = -1;
		int bestpath = -1;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			//2.find vertices that are not yet included in the solution
			int endj = (int)sol.solution[d].size();
			for (int j = 1; j < endj - 1; ++j)// for all  inlcuded vertices in the solution (non-depots)
			{//4. check for interesting replacements
				Ins::Vertex* z = sol.solution[d][j];//z point that will be replaced
				int breakz = sol.action[d][j];
				int nb_size = (int)z->nb[d].size();
				for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
				{
					Ins::Vertex* y = z->nb[d][i];//y replacement
					Ins::Vertex* x = sol.solution[d][j - 1];//predecessor of z
					Ins::Vertex* w = sol.solution[d][j + 1];//successor of z
					int breakw = sol.action[d][j + 1];
					if ((sol.available[y->index] * y->score > z->score) && (x->nbi[d][y->index]) && (y->nbi[d][w->index]))//if interesting and possible
					{
						// check if replacement is feasible
						//gather departure time
						double currenttime = sol.traveltime[d][j - 1] + ins->t[d].EDT;//service bij x zit hier al in
						//travel time from x to y
						double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
						if (arrivaltime + sol.action[d][j] * ins->breakdur < y->LTW[d])
						{
							arrivaltime = y->LTW[d] - sol.action[d][j] * ins->breakdur;
						}
						if (arrivaltime > y->UTW[d])
						{
							continue;//mag je al stoppen met rekenen voor dit punt op deze positie
						}
						arrivaltime += y->serv + breakz * ins->breakdur;
						double ttxy = arrivaltime - time_periods[0];
						//travel time from y to w
						arrivaltime = ins->arrival_time(y->con[w->index], arrivaltime);
						if (arrivaltime + sol.action[d][j + 1] * ins->breakdur < w->LTW[d])
						{
							arrivaltime = w->LTW[d] - sol.action[d][j + 1] * ins->breakdur;
						}
						arrivaltime += w->serv + sol.action[d][j + 1] * ins->breakdur;
						double increase = arrivaltime - currenttime;//new traveltime=>service time included ttxy +ttyw
						increase -= (sol.traveltime[d][j + 1] - sol.traveltime[d][j - 1]);//substract old traveltime=> ttxz + ttzw
						if (increase <= sol.max_shift[d][j + 1])//check of het punt gereplaced kan worden
						{
							improvement = true;
							double ratiocheck;
							if (increase <= 0)
							{
								ratiocheck = y->score - z->score;
							}
							else
							{
								ratiocheck = double(y->score - z->score) / increase;

							}
							if (ratiocheck > ratio)//enkel op minimale increase checken
							{
								ratio = ratiocheck;
								position = j;
								candidate = y;
								ttxybest = ttxy;
								bestpath = d;
							}
						}//end if feasible replacement
					}//if interesting
				}//for i
			}//end for j
		}//end for d
		if (improvement)
		{
			// execute replacement
			sol.replacevertex(bestpath, candidate, position);
			sol.check();
			cout << "debug here" << endl;
		}//end if improvement
	}//end while improvement
}

void Moves::exchange(Sol& sol)
{
	bool improvement = true;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 0;
		int bestpatha;
		int bestpathb;
		int bestindexa;
		int bestindexb;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			for (int i = 0; i < (int)sol.solution[d].size() - 2; ++i)
			{
				Ins::Vertex* x = sol.solution[d][i + 1];
				Ins::Vertex* y = sol.solution[d][i + 2];
				double oldtraveltime = (sol.traveltime[d][i + 2] - (y->serv + sol.action[d][i + 2] * ins->breakdur + x->serv + sol.action[d][i + 1] * ins->breakdur)) - sol.traveltime[d][i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					if (d != e)
					{
						for (int j = 0; j < sol.solution[e].size() - 1; ++j)
						{
							Ins::Vertex* a = sol.solution[e][j];
							Ins::Vertex* b = sol.solution[e][j + 1];
							double departuretime = ins->t[e].EDT + sol.traveltime[e][j];
							double newtraveltime = 0;
							//tw check +nb check
							//traveltime a to x
							double arrivaltime = ins->arrival_time(a->con[x->index], departuretime);
							if (arrivaltime < x->LTW[d])
							{
								arrivaltime = x->LTW[d];
							}
							if (arrivaltime > x->UTW[d])
							{
								continue;
							}
							newtraveltime += arrivaltime - departuretime;
							arrivaltime += x->serv;//don't insert break
							departuretime = arrivaltime;
							//traveltime from x to b
							arrivaltime = ins->arrival_time(x->con[b->index], departuretime);
							if (arrivaltime < b->LTW[d])
							{
								arrivaltime = b->LTW[d];
							}
							if (arrivaltime > b->UTW[d])
							{
								continue;
							}
							newtraveltime += arrivaltime - departuretime;
							double decrease = oldtraveltime - newtraveltime;
							if (decrease > bestdecrease)
							{
								bestdecrease = decrease;
								arrivaltime += b->serv + sol.action[e][j + 1] * ins->breakdur;
								double shiftmaster = (arrivaltime - ins->t[e].EDT) - sol.traveltime[e][j + 1];//increase in travel time when inserting
								if ((shiftmaster <= sol.max_shift[e][j + 1]))//point can be inserted
								{
									improvement = true;
									bestpatha = e;
									bestpathb = d;//slave
									bestindexa = j;
									bestindexb = i + 1;//slave
								}
							}
						}
					}//end for different path
				}//end for slavepath e
			}//end for master sol i
		}// end for masterpath d
		if (improvement)
		{
			//insert new vertex into slave path
			bool insert = false;
			double ratio = 0.0;
			int position = -1;
			Ins::Vertex* x = sol.solution[bestpathb][bestindexb - 1];//point before insertion
			Ins::Vertex* exch = sol.solution[bestpathb][bestindexb];//point that will be exchanged
			sol.scores[bestpathb] -= exch->score;
			sol.volume[bestpathb] -= exch->volume;
			sol.weight[bestpathb] -= exch->weight;
			Ins::Vertex* candidate = NULL;
			int nb_size = (int)x->nb[bestpathb].size();
			for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex
			{
				Ins::Vertex* y = x->nb[bestpathb][i];//point that might be inserted
				Ins::Vertex* z = sol.solution[bestpathb][bestindexb + 1];//point to shift
				int breakz = sol.action[bestpathb][bestindexb + 1];
				if ((sol.available[y->index]) && (y->nbi[bestpathb][z->index]))
				{
					if ((sol.volume[bestpathb] + y->nb[bestpathb][i]->volume < ins->t[bestpathb].V_max) && (sol.weight[bestpathb] + y->nb[bestpathb][i]->weight < ins->t[bestpathb].W_max))
					{
						//gather departure time
						double currenttime = sol.traveltime[bestpathb][bestindexb - 1] + ins->t[bestpathb].EDT;//service bij x zit hier al in
						//travel time from x to y
						double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
						if (arrivaltime < y->LTW[bestpathb])
						{
							arrivaltime = y->LTW[bestpathb];
						}
						if (arrivaltime > y->UTW[bestpathb])
						{
							continue;//mag je al stoppen met rekenen voor dit punt op deze positie
						}
						arrivaltime += y->serv + sol.action[bestpathb][bestindexb] * ins->breakdur;
						//travel time from y to z
						arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
						if (arrivaltime < z->LTW[bestpathb])
						{
							arrivaltime = z->LTW[bestpathb];
						}
						arrivaltime += z->serv + breakz * ins->breakdur;
						double shift = (arrivaltime - ins->t[bestpathb].EDT) - sol.traveltime[bestpathb][bestindexb + 1];//increase in travel time
						if (shift <= sol.max_shift[bestpathb][bestindexb + 1])//check of het punt geinsert kan worden
						{
							double ratiocheck;
							if (shift <= 0)
							{
								ratiocheck = y->score;
							}
							else
							{
								//ratiocheck = 1.0 / shift;
								ratiocheck = double(y->score * y->score) / shift;
								//score/shift in traveltime
							}
							if (ratiocheck > ratio)//enkel op minimale increase checken
							{
								insert = true;
								ratio = ratiocheck;
								position = bestindexb;
								candidate = y;
							}
							//update candidates
						}// end if check feasible insertion
					}//if cap restrictions
				}//end if available
			}//end i
			if (insert)
			{
				// execute best insertion if improvement found
				sol.available[candidate->index] = false;
				sol.solution[bestpathb][bestindexb] = candidate;
				sol.scores[bestpathb] += candidate->score;// update score of the new solution
				sol.score += candidate->score;// update score of the new solution
				sol.volume[bestpathb] += candidate->volume;
				sol.weight[bestpathb] += candidate->weight;
				//update travel time solution vector 
				double currenttime = sol.traveltime[bestpathb][bestindexb - 1] + ins->t[bestpathb].EDT;
				for (int u = position; u < sol.solution[bestpathb].size(); ++u)
				{
					//gather departure time and corresponding time slot
					Ins::Vertex* o = sol.solution[bestpathb][u - 1];
					Ins::Vertex* p = sol.solution[bestpathb][u];
					//travel time from van o to p
					double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
					if (arrivaltime < p->LTW[bestpathb])
					{
						arrivaltime = p->LTW[bestpathb];
					}
					arrivaltime += p->serv + (sol.action[bestpathb][u] * ins->breakdur);
					sol.max_shift[bestpathb][u] = (sol.traveltime[bestpathb][u] + sol.max_shift[bestpathb][u]) - (arrivaltime - ins->t[bestpathb].EDT);
					//old arrival time + old maxshift - new arrival time = new max_shift
					sol.traveltime[bestpathb][u] = arrivaltime - ins->t[bestpathb].EDT;
					currenttime = arrivaltime;
				}//end for
				//the value of max_shift j+1 is now wrong but will be soon be updated
				//complete update of max_shift for the included vertices before the insertion
				double departuretime = 0;
				double arrivaltime = (sol.traveltime[bestpathb][position + 1] + ins->t[bestpathb].EDT + sol.max_shift[bestpathb][position + 1]) - sol.solution[bestpathb][position + 1]->serv + sol.action[bestpathb][position + 1] * ins->breakdur;//service time eraftrekken
				for (int vz = position; vz > 0; --vz)
				{
					//define the 2 elements
					Ins::Vertex* f = sol.solution[bestpathb][vz];
					int breakf = sol.action[bestpathb][vz];
					Ins::Vertex* g = sol.solution[bestpathb][vz + 1];
					//find proper departure time and time slot
					departuretime = ins->departure_time(f->con[g->index], arrivaltime);	//utw check
					if (departuretime > f->UTW[bestpathb] + f->serv + breakf * ins->breakdur)
					{
						departuretime = f->UTW[bestpathb] + f->serv + breakf * ins->breakdur;
					}
					if (breakf == 1)
					{
						if (departuretime > ins->breakend + ins->breakdur)
						{
							departuretime = ins->breakend + ins->breakdur;
						}
					}
					sol.max_shift[bestpathb][vz] = departuretime - (sol.traveltime[bestpathb][vz] + ins->t[bestpathb].EDT);
					departuretime -= f->serv + breakf * ins->breakdur;
					arrivaltime = departuretime;
				}// end for vz

				//update master
				sol.solution[bestpatha].insert(sol.solution[bestpatha].begin() + bestindexa + 1, exch);//insert exch after x
				sol.traveltime[bestpatha].insert(sol.traveltime[bestpatha].begin() + bestindexa + 1, 0);//insert temporary value
				sol.action[bestpatha].insert(sol.action[bestpatha].begin() + bestindexa + 1, 0);//insert regular visit action change later when necessary
				if (bestindexa < sol.breakindex[bestpatha])
				{
					sol.breakindex[bestpatha] += 1;//update breakindex
				}
				sol.max_shift[bestpatha].insert(sol.max_shift[bestpatha].begin() + bestindexa + 1, 0);
				sol.scores[bestpatha] += exch->score;// update score of the new solution
				//score of solution remains the same
				sol.volume[bestpatha] += exch->volume;
				sol.weight[bestpatha] += exch->weight;
				//update travel time solution vector 
				currenttime = sol.traveltime[bestpatha][bestindexa] + ins->t[bestpatha].EDT;
				for (int u = bestindexa + 1; u < sol.solution[bestpatha].size(); ++u)
				{
					//gather departure time and corresponding time slot
					Ins::Vertex* o = sol.solution[bestpatha][u - 1];
					Ins::Vertex* p = sol.solution[bestpatha][u];
					//travel time from van o to p
					double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
					if (arrivaltime < p->LTW[bestpatha])
					{
						arrivaltime = p->LTW[bestpatha];
					}
					arrivaltime += p->serv + sol.action[bestpatha][u] * ins->breakdur;
					sol.max_shift[bestpatha][u] = (sol.traveltime[bestpatha][u] + sol.max_shift[bestpatha][u]) - (arrivaltime - ins->t[bestpatha].EDT);
					//old arrival time + old maxshift - new arrival time = new max_shift
					sol.traveltime[bestpatha][u] = arrivaltime - ins->t[bestpatha].EDT;
					currenttime = arrivaltime;
				}//end for
				//the value of max_shift j+1 is now wrong but will be soon be updated
				//complete update of max_shift for the included vertices before the insertion
				departuretime = 0;
				arrivaltime = (sol.traveltime[bestpatha][bestindexa + 2] + ins->t[bestpatha].EDT + sol.max_shift[bestpatha][bestindexa + 2]) - sol.solution[bestpatha][bestindexa + 2]->serv + sol.action[bestpatha][bestindexa + 2] * ins->breakdur;//service time eraftrekken
				for (int vz = bestindexa + 1; vz > 0; --vz)
				{
					//define the 2 elements
					Ins::Vertex* f = sol.solution[bestpatha][vz];
					int breakf = sol.action[bestpatha][vz];
					Ins::Vertex* g = sol.solution[bestpatha][vz + 1];
					//find proper departure time and time slot
					departuretime = ins->departure_time(f->con[g->index], arrivaltime);	//utw check
					if (departuretime > f->UTW[bestpatha] + f->serv + breakf * ins->breakdur)
					{
						departuretime = f->UTW[bestpatha] + f->serv + breakf * ins->breakdur;
					}
					if (breakf == 1)
					{
						if (departuretime > ins->breakend + ins->breakdur)
						{
							departuretime = ins->breakend + ins->breakdur;
						}
					}
					sol.max_shift[bestpatha][vz] = departuretime - (sol.traveltime[bestpatha][vz] + ins->t[bestpatha].EDT);
					departuretime -= f->serv + breakf * ins->breakdur;
					arrivaltime = departuretime;
				}// end for vz
			}
			else
			{
				improvement = false;
				sol.scores[bestpathb] += exch->score;
				sol.volume[bestpathb] += exch->volume;
				sol.weight[bestpathb] += exch->weight;
			}
			sol.check();
			cout << "hier" << endl;
		}//end if improvement
	}//end while improvement

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
							arrivaltime += second->serv + sol.action[d][l + 1] * ins->breakdur;
							currenttime = arrivaltime;
						}// end l
						newtraveltime = currenttime - departuretime;
						oldtraveltime = ((ins->t[d].EDT + sol.traveltime[d][j + 1]) - departuretime);
						delta_tt = newtraveltime - oldtraveltime;
						if (delta_tt < 0)
						{
							double remember = sol.traveltime[d].back();
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
								cout << "actual improvement" << remember - sol.traveltime[d].back() << endl;
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

void Moves::two_opt_nb(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		bool improvement = true;
		while (improvement)
		{
			improvement = false;
			int end = (int)sol.solution[d].size();
			for (int i = 1; i < end - 1; ++i)// depots don't count
			{
				for (int j = i + 1; j < end - 2; ++j)
				{
					Ins::Vertex* k = sol.solution[d][i - 1];//pred partner 1
					Ins::Vertex* y = sol.solution[d][i];//partner 1
					Ins::Vertex* l = sol.solution[d][i + 1];//suc partner 1
					Ins::Vertex* m = sol.solution[d][j - 1];//pred  partner 2
					Ins::Vertex* z = sol.solution[d][j];// partner 2
					Ins::Vertex* n = sol.solution[d][j + 1];//suc partner 2
					if ((k->nbi[d][z->index]) && (z->nbi[d][m->index]) && (l->nbi[d][y->index]) && (y->nbi[d][n->index]))
					{
						//tijdelijk opten
						for (int f = 0; f < 1 + (j - i) / 2; ++f)
						{
							Ins::Vertex* temp = sol.solution[d][j - f];
							sol.solution[d][j - f] = sol.solution[d][i + f];
							sol.solution[d][i + f] = temp;
						}
						bool reval = false;
						double delta_tt = 0;
						double departuretime = ins->t[d].EDT + sol.traveltime[d][i - 1];//bij punt voor y
						double currenttime = departuretime;
						double arrivaltime;
						for (int l = i - 1; l < j + 1; ++l)
						{
							Ins::Vertex* first = sol.solution[d][l];
							Ins::Vertex* second = sol.solution[d][l + 1];
							arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
							if (sol.action[d][l + 1] == 1)
							{
								if ((arrivaltime < ins->breakstart) || (arrivaltime > ins->breakend))
								{//break is wrongly positioned in local evaluation
									reval = true;
								}
							}
							if (arrivaltime + sol.action[d][l + 1] * ins->breakdur < second->LTW[d])
							{
								arrivaltime = second->LTW[d] - sol.action[d][l + 1] * ins->breakdur;
							}
							if (arrivaltime > second->UTW[d])
							{
								goto stop;//swap is infeasible, stop for this point
							}
							//als hij boven de utw gaat dan is delta_tt toch negatief
							currenttime = arrivaltime + second->serv + sol.action[d][l + 1] * ins->breakdur;
						}// end l
						//int remember = sol.traveltime[d].back();
						delta_tt = currenttime - (ins->t[d].EDT + sol.traveltime[d][j + 1]);
						//break wijzigingen kunnen niet bij een locale evaluatie
						if (delta_tt < 0)//local feasible 2-opt gevonden 
						{
							//let op delta is een negatief getal
							if ((sol.breakindex[d] > j + 1) && ((ins->t[d].EDT + sol.traveltime[d][sol.breakindex[d]]) - ((sol.solution[d][sol.breakindex[d]]->serv + ins->breakdur) - delta_tt) < ins->breakstart))
							{
								reval = true;
							}

							if (reval)
							{//global evalution necessary if the breaks need to be repositioned
								bool reqbreak = true;
								vector<double> temptravel;
								temptravel.push_back(0);
								int breakindex = -1;
								currenttime = ins->t[d].EDT;
								for (int vv = 0; vv < end - 1; ++vv)
								{
									Ins::Vertex* first = sol.solution[d][vv];
									Ins::Vertex* second = sol.solution[d][vv + 1];
									int breaksecond = 0;
									arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
									if ((reqbreak) && (arrivaltime >= ins->breakstart))
									{
										breaksecond = 1;
										breakindex = vv + 1;
										reqbreak = false;
									}
									if (arrivaltime + breaksecond * ins->breakdur < second->LTW[d])
									{
										arrivaltime = second->LTW[d] - breaksecond * ins->breakdur;
									}
									arrivaltime += second->serv + breaksecond * ins->breakdur;
									temptravel.push_back(arrivaltime - ins->t[d].EDT);
									currenttime = arrivaltime;
								}
								if (reqbreak == true)//als je aankomt bij het einddepot en nog steeds geen break genomen hebt
								{
									breakindex = end - 1;
									currenttime += ins->breakdur;//breaktime bijtellen bij aankomst tijd bij einddepot
								}
								if (currenttime - ins->t[d].EDT >= sol.traveltime[d].back())
								{//solutiond does not improve
									goto stop;
								}
								else
								{//solution improves
									improvement = true;
									sol.action[d][sol.breakindex[d]] = 0;
									sol.action[d][breakindex] = 1;
									sol.breakindex[d] = breakindex;
									sol.traveltime[d] = temptravel;
									break;
								}
							}
							else
							{//geen global reval nodig want break blijf op dezelfde plaats
								improvement = true;
								//travel time herevalueren
								currenttime = departuretime;
								for (int tt = i; tt < end; ++tt)
								{
									Ins::Vertex* first = sol.solution[d][tt - 1];
									Ins::Vertex* second = sol.solution[d][tt];
									int breakz = sol.action[d][tt];
									arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
									//add possible waiting time
									if (arrivaltime + sol.action[d][tt] * ins->breakdur < second->LTW[d])
									{
										arrivaltime = second->LTW[d] - sol.action[d][tt] * ins->breakdur;
									}
									arrivaltime += second->serv + breakz * ins->breakdur;
									//sol.max_shift[d][tt] = (sol.traveltime[d][tt] + sol.max_shift[d][tt]) - (arrivaltime - t[d].EDT);
									sol.traveltime[d][tt] = arrivaltime - ins->t[d].EDT;
									currenttime = arrivaltime;
								}//end for t
								//cout<<"actual improvement"<<sol.traveltime[d].back()-remember<<endl;
								//revaluate max_shift for all vertices before j+1
								/*
								int arrivaltime = sol.traveltime[d][j + 1] + t[d].EDT + sol.max_shift[d][j + 1] - sol.solution[d][j + 1]->serv[sol.action[d][j+1]];
								for (int tt = j; tt > 0; --tt)
								{
								//define the 2 elements
								Vertex *first = sol.solution[d][tt];
								int breakz = sol.action[d][tt];
								Vertex *second = sol.solution[d][tt + 1];
								//find proper departure time
								departuretime = departure_time(first->con[second->index], arrivaltime);
								if (departuretime > first->UTW[day] + first->serv[breakz])//utw check
								{
								departuretime = first->UTW[day] + first->serv[breakz];
								}
								//store result
								sol.max_shift[d][tt] = departuretime - (sol.traveltime[d][tt] + t[d].EDT);
								//update arrivaltime for the calculation of next point
								arrivaltime = departuretime - first->serv[breakz];
								}// end for t
								*/
								sol.check();
								cout << "hier" << endl;
								break;
								//i is geswapped ofwel j lus terug doorlopen ofwel naar volgende i gaan
							}//end if feas

						}//end executed 2 opt
						else
						{
						stop:
							//terug zetten
							for (int f = 0; f < 1 + (j - i) / 2; ++f)
							{
								Ins::Vertex* temp = sol.solution[d][j - f];
								sol.solution[d][j - f] = sol.solution[d][i + f];
								sol.solution[d][i + f] = temp;
							}
						}
					}
				}// end for al j
			}// end for all i
		}
	}//end for all paths
}

void Moves::best_move_nb(Sol& sol)//move 1 vertex from one path to another to save travel time
{
	bool improvement = true;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 0;
		int bestd;
		int beste;
		int besti;
		int bestj;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			for (int i = 0; i < sol.solution[d].size() - 2; ++i)
			{
				Ins::Vertex* w = sol.solution[d][i];
				Ins::Vertex* x = sol.solution[d][i + 1];
				Ins::Vertex* y = sol.solution[d][i + 2];
				int breaky = sol.action[d][i + 2];
				double ttwxy = (sol.traveltime[d][i + 2] - (y->serv + sol.action[d][i + 2] * ins->breakdur + x->serv + sol.action[d][i + 1] * ins->breakdur)) - sol.traveltime[d][i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					if (d != e)//no relocate on the same path
					{
						for (int j = 0; j < sol.solution[e].size() - 2; ++j)
						{
							Ins::Vertex* a = sol.solution[e][j];
							Ins::Vertex* b = sol.solution[e][j + 1];
							if ((a->nbi[e][x->index]) && (x->nbi[e][b->index]) && (w->nbi[d][y->index]))
							{
								if ((sol.weight[e] + x->weight < ins->t[e].W_max) && (sol.volume[e] + x->volume < ins->t[e].V_max))
								{
									bool reqbreak = false;
									if (sol.breakindex[d] >= i + 1)//als de break op x staat moet hij verschoven worden naar achter+ door het onstane gat kan de break ook te vroeg komen en dan moet hij ook naar achter verschoven worden
									{
										reqbreak = true;
									}
									//local evaluation on path e
									double ttab = (sol.traveltime[e][j + 1] - (b->serv + sol.action[e][j + 1] * ins->breakdur)) - sol.traveltime[e][j];
									//calculate axb
									double departuretime = ins->t[e].EDT + sol.traveltime[e][j];
									//traveltime a to x
									double arrivaltime = ins->arrival_time(a->con[x->index], departuretime);
									if (arrivaltime < x->LTW[e])//je mag niet breaken op x want er is al een break op de e route
									{
										arrivaltime = x->LTW[e];
									}
									if (arrivaltime > x->UTW[e])
									{
										continue;
									}
									double ttaxb = arrivaltime - departuretime;
									//insert vertex without break
									arrivaltime += x->serv;
									departuretime = arrivaltime;
									//traveltime from x to b
									arrivaltime = ins->arrival_time(x->con[b->index], departuretime);
									if (arrivaltime + sol.action[e][j + 1] * ins->breakdur < b->LTW[e])
									{
										arrivaltime = b->LTW[e] - (sol.action[e][j + 1] * ins->breakdur);
									}
									if (arrivaltime > b->UTW[e])
									{
										continue;
									}
									ttaxb += arrivaltime - departuretime;
									double arrivalb = arrivaltime + b->serv + sol.action[e][j + 1] * ins->breakdur;
									//on path d now
									//calculate wy
									departuretime = ins->t[d].EDT + sol.traveltime[d][i];
									//traveltime w to y
									arrivaltime = ins->arrival_time(w->con[y->index], departuretime);
									//break niet meetellen voor local evaluation
									if ((reqbreak) && (arrivaltime >= ins->breakstart))
									{//break na y wordt bij global evaluation in rekening gebracht
										reqbreak = false;
										breaky = 1;
										//de oude 1 op x gaat verwijderd worden
									}
									if (arrivaltime + breaky * ins->breakdur < y->LTW[d])
									{
										arrivaltime = y->LTW[d] - (breaky * ins->breakdur);
									}
									if (arrivaltime > y->UTW[d])
									{
										continue;
									}
									double ttwy = arrivaltime - departuretime;
									double arrivaly = arrivaltime + y->serv + breaky * ins->breakdur;
									double localincrease = (ttaxb - ttab) + x->serv;
									double localdecreasetotal = ttwxy + ttab - (ttwy + ttaxb);
									//local improvement check
									if ((localincrease <= sol.max_shift[e][j + 1]) && (localdecreasetotal > bestdecrease))
									{
										//global improvement check
										//check enddepot time on path d
										double currenttime = arrivaly;
										for (int m = i + 3; m < sol.solution[d].size(); ++m)
										{
											//gather departure time and corresponding time slot
											Ins::Vertex* o = sol.solution[d][m - 1];
											Ins::Vertex* p = sol.solution[d][m];
											int breakp = sol.action[d][m];
											//travel time from van o to p
											double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
											if ((reqbreak) && (arrivaltime >= ins->breakstart))
											{
												breakp = 1;
												reqbreak = false;
											}
											if (arrivaltime + breakp * ins->breakdur < p->LTW[d])
											{
												arrivaltime = p->LTW[d] - (breakp * ins->breakdur);
											}
											arrivaltime += p->serv + breakp * ins->breakdur;
											currenttime = arrivaltime;
										}
										if (reqbreak)//nog altijd geen break kunnen plaatsen dan moet hij op het einddepot komen
										{
											currenttime += ins->breakdur;
										}
										double globaldecreasetotal = (sol.traveltime[d].back() - (currenttime - ins->t[d].EDT));
										//cout << "enddepot time d: " << currenttime - t[d].EDT << endl;
										//cout << "decrease d: " << (sol.traveltime[d].back() - (currenttime - t[d].EDT)) << endl;
										//check enddepot time on path e
										currenttime = arrivalb;
										for (int m = j + 2; m < (int)sol.solution[e].size(); ++m)
										{
											//gather departure time and corresponding time slot
											Ins::Vertex* o = sol.solution[e][m - 1];
											Ins::Vertex* p = sol.solution[e][m];
											//travel time from van o to p
											double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
											if (arrivaltime + sol.action[e][m] * ins->breakdur < p->LTW[e])
											{
												arrivaltime = p->LTW[e] - (sol.action[e][m] * ins->breakdur);
											}
											arrivaltime += p->serv + sol.action[e][m] * ins->breakdur;
											currenttime = arrivaltime;
										}
										//cout << "enddepot time e: " << currenttime - t[e].EDT << endl;
										//cout << "increase e: " << (currenttime - t[e].EDT) - sol.traveltime[e].back() << endl;
										globaldecreasetotal -= ((currenttime - ins->t[e].EDT) - sol.traveltime[e].back());
										if (globaldecreasetotal > bestdecrease)
										{
											bestdecrease = globaldecreasetotal;
											beste = e;
											bestj = j;
											bestd = d;
											besti = i + 1;
											improvement = true;
										}
									}
								}//end if still capacity free
							}//end if neighbour
						}//end for j
					}//end for d!=e
				}//end for path e
			}//end for master sol i
		}// end for masterpath d
		if (improvement)
		{
			Ins::Vertex* candidate = sol.solution[bestd][besti];
			sol.removevertex(bestd, besti);//remove vertex from path d
			if ((sol.breakindex[bestd] == int(sol.solution.size()) - 1) && (ins->t[bestd].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "break pulled" << endl;
				pull_break(sol, bestd);
			}
			sol.insertvertex(beste, candidate, bestj);//insert vertex on path e
			if ((sol.breakindex[beste] == int (sol.solution[beste].size()) - 1) && (ins->t[beste].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "break pulled" << endl;
				pull_break(sol, beste);
			}
		}//end if improvement
	}//end while improvement
}

void Moves::best_swap2_nb(Sol& sol)
{
	bool improvement = true;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 0;
		int bestd;
		int beste;
		int besti;
		int bestj;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			for (int i = 0; i < sol.solution[d].size() - 2; ++i)
			{
				Ins::Vertex* w = sol.solution[d][i];
				Ins::Vertex* x = sol.solution[d][i + 1];
				Ins::Vertex* y = sol.solution[d][i + 2];
				int breaky = sol.action[d][i + 2];
				double ttwxy = (sol.traveltime[d][i + 2] - (y->serv + sol.action[d][i + 2] * ins->breakdur + x->serv + sol.action[d][i + 1] * ins->breakdur)) - sol.traveltime[d][i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					if (d != e)//no relocate on the same path
					{
						for (int j = 0; j < sol.solution[e].size() - 2; ++j)
						{
							Ins::Vertex* a = sol.solution[e][j];
							Ins::Vertex* b = sol.solution[e][j + 1];
							Ins::Vertex* c = sol.solution[e][j + 2];
							if ((a->nbi[e][x->index]) && (x->nbi[e][b->index]) && (w->nbi[d][y->index]))
							{
								if ((sol.weight[e] + (x->weight-b->weight) < ins->t[e].W_max) && (sol.volume[e] + (x->weight-b->volume) < ins->t[e].V_max)&&(sol.weight[d] + (b->weight - x->weight) < ins->t[d].W_max)&&(sol.volume[d] + (b->volume - x->volume) < ins->t[d].V_max))
								{//check if a potential increase in volume and weight is allowed on each tour
									
									//local evaluation on path d
									double departuretime = ins->t[d].EDT + sol.traveltime[d][i];
									// w to b
									double arrivaltime = ins->arrival_time(w->con[b->index], departuretime);
									if (arrivaltime + (sol.action[d][i + 1] * ins->breakdur) < b->LTW[d])
									{
										arrivaltime = b->LTW[d];
									}
									if (arrivaltime > b->UTW[d])
									{
										continue;
									}
									arrivaltime += b->serv+(sol.action[d][i+1]*ins->breakdur);
									departuretime = arrivaltime;
									//traveltime from b to y
									arrivaltime = ins->arrival_time(b->con[y->index], departuretime);
									if (arrivaltime + sol.action[d][i + 2] * ins->breakdur < b->LTW[d])
									{
										arrivaltime = b->LTW[d] - (sol.action[d][i + 2] * ins->breakdur);
									}
									if (arrivaltime > b->UTW[d])
									{
										continue;
									}
									double arrivaltimey =arrivaltime+ y->serv + (sol.action[d][i + 2] * ins->breakdur);
									double diffd = (ins->t[d].EDT+sol.traveltime[d][i + 2])-arrivaltimey;

									//local evaluation on path e
									departuretime = ins->t[e].EDT + sol.traveltime[e][j];
									//traveltime a to x
									arrivaltime = ins->arrival_time(a->con[x->index], departuretime);
									if (arrivaltime < x->LTW[e])//je mag niet breaken op x want er is al een break op de e route
									{
										arrivaltime = x->LTW[e];
									}
									if (arrivaltime > x->UTW[e])
									{
										continue;
									}
									arrivaltime += x->serv + (sol.action[e][j + 1] * ins->breakdur);
									departuretime = arrivaltime;
									//traveltime from x to c
									arrivaltime = ins->arrival_time(x->con[c->index], departuretime);
									if (arrivaltime + sol.action[e][j + 2] * ins->breakdur < c->LTW[e])
									{
										arrivaltime = c->LTW[e] - (sol.action[e][j + 2] * ins->breakdur);
									}
									if (arrivaltime > c->UTW[e])
									{
										continue;
									}
									double arrivaltimec =arrivaltime+ c->serv + (sol.action[e][j + 2] * ins->breakdur);
									double diffe = (ins->t[e].EDT+sol.traveltime[e][j + 2]) - arrivaltimec;
									double localdecreasetotal = diffd + diffe;
									//local improvement check: check if potential increase is allowed and whether there is a overall travel time gain
									if ((diffd<= sol.max_shift[d][i + 1]) && (diffe <= sol.max_shift[e][j + 1]) && (localdecreasetotal > bestdecrease))
									{
										//global improvement check
										//check enddepot time on path d
										double currenttime = arrivaltimey;
										for (int m = i + 3; m < (int) sol.solution[d].size(); ++m)
										{
											//gather departure time and corresponding time slot
											Ins::Vertex* o = sol.solution[d][m - 1];
											Ins::Vertex* p = sol.solution[d][m];
											//travel time from van o to p
											double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
											if (arrivaltime + (sol.action[d][m] * ins->breakdur) < p->LTW[d])
											{
												arrivaltime = p->LTW[d] - (sol.action[d][m] * ins->breakdur);
											}
											arrivaltime += p->serv + (sol.action[d][m] * ins->breakdur);
											currenttime = arrivaltime;
										}
										double globaldecreasetotal = (sol.traveltime[d].back() - (currenttime - ins->t[d].EDT));
										//cout << "enddepot time d: " << currenttime - t[d].EDT << endl;
										//cout << "decrease d: " << (sol.traveltime[d].back() - (currenttime - t[d].EDT)) << endl;
										//check enddepot time on path e
										currenttime = arrivaltimec;
										for (int m = j + 3; m < (int)sol.solution[e].size(); ++m)
										{
											//gather departure time and corresponding time slot
											Ins::Vertex* o = sol.solution[e][m - 1];
											Ins::Vertex* p = sol.solution[e][m];
											//travel time from van o to p
											double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
											if (arrivaltime + sol.action[e][m] * ins->breakdur < p->LTW[e])
											{
												arrivaltime = p->LTW[e] - (sol.action[e][m] * ins->breakdur);
											}
											arrivaltime += p->serv + sol.action[e][m] * ins->breakdur;
											currenttime = arrivaltime;
										}
										//cout << "enddepot time e: " << currenttime - t[e].EDT << endl;
										//cout << "increase e: " << (currenttime - t[e].EDT) - sol.traveltime[e].back() << endl;
										globaldecreasetotal -= ((currenttime - ins->t[e].EDT) - sol.traveltime[e].back());
										if (globaldecreasetotal > bestdecrease)
										{
											bestdecrease = globaldecreasetotal;
											beste = e;
											bestj = j+1;
											bestd = d;
											besti = i + 1;
											improvement = true;
										}
									}
								}//end if still capacity free
							}//end if neighbour
						}//end for j
					}//end for d!=e
				}//end for path e
			}//end for master sol i
		}// end for masterpath d
		if (improvement)
		{
			Ins::Vertex *x = sol.solution[bestd][besti];
			Ins::Vertex *b = sol.solution[beste][bestj];
			sol.replacevertex(bestd, b, besti);
			sol.replacevertex(beste, x, bestj);
		}//end if improvement
	}//end while improvement

}
