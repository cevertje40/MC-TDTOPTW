#include "moves.h"

bool Moves::insert_nb(Sol& sol, int mode)//insert vertex into a tour in order to increase the score
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestratio = 0.0;
		int position = -1;
		Ins::Vertex* candidate = NULL;
		Sol::Tour* besttour = NULL;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			int t = sol.tourindex[d];
			Sol::Tour* tour= &sol.tours[sol.tourindex[d]];
			//sort possible candidates
			int endj = (int)sol.tours[t].seq.size();
			for (int j = 0; j < endj - 1; ++j)// for all  inlcuded vertices in the solution (non-depot)
			{
				Ins::Vertex* x = tour->seq[j];//point before insertion
				int nb_size = (int)x->nb[t].size();
				for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex
				{
					Ins::Vertex* y = x->nb[t][i];//point that might be inserted
					Ins::Vertex* z = tour->seq[j + 1];//point to shift
					int breakz = tour->action[j + 1];
					if ((sol.available[y->index]) && (y->nbi[t][z->index]))//y moet buur van z zijn want 
					{
						if ((tour->volume + y->volume <= ins->t[t].V_max) && (tour->weight + y->weight <= ins->t[t].W_max))//check capacity constraint
						{
							//gather departure time
							double currenttime = tour->deptime[j] + ins->t[t].EDT;//service bij x zit hier al in
							//travel time from x to y
							double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
							if (arrivaltime < y->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
							{
								arrivaltime = y->LTW[t];
							}
							if (arrivaltime > y->UTW[t])
							{
								continue;//mag je al stoppen met rekenen voor dit punt op deze positie
							}
							arrivaltime += y->serv;//break inserten kan niet
							//travel time from y to z
							double waitz = 0;//wait at z
							arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
							if (arrivaltime + breakz * (ins->breakdur) < z->LTW[t])
							{
								waitz = z->LTW[t] - (arrivaltime + breakz * ins->breakdur);
								arrivaltime = z->LTW[t] - (breakz * ins->breakdur);

							}
							arrivaltime += z->serv + breakz * ins->breakdur;
							double shift = (arrivaltime - ins->t[t].EDT) - tour->deptime[j + 1];//increase in travel time
							if (shift <= tour->max_shift[j + 1])//check of het punt geinsert kan worden
							{
								
								//double weightavail = ins->t[t].W_max - tour->weight;
								//double volumeavail = ins->t[t].V_max - tour->volume;
								//double ratiocheck;
								//double consumption = 1;
								//if (shift <= waitz)//shift smaller than wait or 0 means there is no time consumption
								//{//avg consumption of 2 resources
									//consumption = ((y->weight / weightavail) + (y->volume / volumeavail))/2;
								//}
								//else
								//{//avg consumption of 3 resources
									//consumption = ((shift / tour->max_shift[j + 1]) + (y->weight / weightavail) + (y->volume / volumeavail))/3;
								//}
								//ratiocheck = (double(y->score * y->score) / consumption);
								double ratiocheck = y->score / y->weight;
								if (ratiocheck > bestratio)//enkel op minimale increase checken
								{//update candidates
									improvement = true;
									bestratio = ratiocheck;
									position = j;
									candidate = y;
									besttour = tour;
									if (mode == 0)//first improvement, otherwise best improvement
									{
										goto insert;
									}
								}
							}// end if check feasible insertion
						}//if cap restrictions
					}//end if available
				}//end i
			}// end for j
		}//end for al paths
		if (improvement)
		{
		insert:
			Sol remember = sol;
			sol.insert_vertex(*besttour, candidate, position);
			//try to pull break if the break is still positioned at the end depot
			int end = (int) besttour->seq.size();
			if ((besttour->breakindex == end - 1) && (ins->t[besttour->index].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "pull break" << endl;
				pull_break(sol, besttour->index);
			}
			//sol.check();
			succes = true;
		}//end improvement
	}//end while improvement
	return succes;
}

bool Moves::exchange_nb(Sol& sol, int mode)//replace a vertex of a tour with non included vertex in the same position in order to increase the score
{
	//1.repeat until no improvement can be found
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestratio = 0.0;
		int position = -1;
		Ins::Vertex* candidate = NULL;
		double ttxybest = -1;
		Sol::Tour* besttour = NULL;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			int t = sol.tourindex[d];//pick a random tour
			Sol::Tour* tour = &sol.tours[sol.tourindex[d]];
			//find vertices that are not yet included in the solution
			int endj = (int)tour->seq.size();
			for (int j = 1; j < endj - 1; ++j)// for all  inlcuded vertices in the solution (non-depots)
			{// check for interesting replacements
				Ins::Vertex* x = tour->seq[j - 1];//predecessor of z
				Ins::Vertex* z = tour->seq[j];//z point that will be replaced
				Ins::Vertex* w = tour->seq[j + 1];//successor of z
				int nb_size = (int)z->nb[t].size();
				for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
				{
					Ins::Vertex* y = z->nb[t][i];//y replacement
					if ((sol.available[y->index] * y->score > z->score) && (x->nbi[t][y->index]) && (y->nbi[t][w->index]))//if interesting and possible
					{
						//add weight and volume check
						if((tour->weight + (y->weight - z->weight) <= ins->t[t].W_max) && (tour->volume + (y->volume - z->volume) <= ins->t[t].V_max))
						{
						
							int breaky = tour->action[j];
							int breakw = tour->action[j+1];
							//gather departure time
							double currenttime = tour->deptime[j - 1] + ins->t[t].EDT;//service bij x zit hier al in
							//travel time from x to y
							double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
							if (arrivaltime +(breaky*ins->breakdur) < y->LTW[t])
							{
								arrivaltime = y->LTW[t]- (breaky*ins->breakdur);
							}
							//break can be scheduled too early or too late
							if (breaky)
							{
								if ((arrivaltime < ins->breakstart)||(arrivaltime>ins->breakend))
									continue;
								else
									arrivaltime += ins->breakdur;
							}
							if (arrivaltime > y->UTW[t])
							{
								continue;//stop the calculation
							}
							arrivaltime += y->serv;
							double ttxy = arrivaltime - ins->t[tour->index].EDT;
							//travel time from y to w
							arrivaltime = ins->arrival_time(y->con[w->index], arrivaltime);
							if (arrivaltime +(breakw*ins->breakdur) < w->LTW[t])
							{
								arrivaltime = w->LTW[t]-(breakw*ins->breakdur);
							}
							//break can be scheduled too early (max_shift check too late)
							if (breakw)
							{
								if (arrivaltime < ins->breakstart)
									continue;
								else
									arrivaltime += ins->breakdur;
							}
							if (arrivaltime > w->UTW[t])
							{
								continue;//stop the calculation
							}
							arrivaltime += w->serv;
							double increase = arrivaltime - currenttime;//new traveltime=>service time included ttxy +ttyw
							increase -= (tour->deptime[j + 1] - tour->deptime[j - 1]);//substract old traveltime=> ttxz + ttzw
							//increase can be negative which could lead to a break that comes to early
							if (increase <= tour->max_shift[j + 1])//check of het punt gereplaced kan worden
							{
								bool breakcheck = true;
								if (increase<0)
								{//global evaluation is necessary as break can fall before breakstart when increase is negative
									double currenttime = arrivaltime;
									for (int m = j+2; m < (int)tour->seq.size(); ++m)
									{
										Ins::Vertex* o = tour->seq[m - 1];
										Ins::Vertex* p = tour->seq[m];//can be the end depot
										int breakp = tour->action[m];
										//travel time from o to p
										double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
										if (arrivaltime + (breakp * ins->breakdur) < p->LTW[t])
										{
											arrivaltime = p->LTW[t] - (breakp * ins->breakdur);
										}
										if (breakp)
										{
											if (arrivaltime < ins->breakstart)
											{
												breakcheck = false;
												continue;
											}
										}
										arrivaltime += p->serv + (breakp * ins->breakdur);
										currenttime = arrivaltime;
									}
								}
								if(breakcheck)
								{//
									
									//double ratiocheck;
									//double weightincrease= (y->weight - z->weight);
									//double weightavailable=ins->t[t].W_max-tour->weight;
									//double volumeincrease= (y->volume - z->volume);
									//double volumeavailable=ins->t[t].V_max-tour->volume;
									//double consumption = 1;
									//if (increase <= 0)//shift smaller than wait or 0 means there is no time consumption
									//{//avg consumption of 2 resources
									//	consumption = (weightincrease/weightavailable) + (volumeincrease/volumeavailable) / 2;
									//}
									//else
									//{//avg consumption of 3 resources
									//	consumption = ((increase/ tour->max_shift[j + 1]) + (weightincrease/weightavailable) + (volumeincrease/volumeavailable)) / 3;
									//}
									//ratiocheck = double(y->score - z->score) / consumption;
									//double ratiocheck = double(y->score - z->score);
									double ratiocheck = double(y->score - z->score) / max(1, (y->weight - z->weight));
									if (ratiocheck > bestratio)//enkel op minimale increase checken
									{
										improvement = true;
										bestratio = ratiocheck;
										position = j;
										candidate = y;
										ttxybest = ttxy;
										besttour = tour;
										if (mode == 0)//first improvement, otherwise best improvement
										{
											goto replace;
										}
									}
								}
							}//end if feasible replacement
						}//capacity checks
					}//if interesting
				}//for i
			}//end for j
		}//end for d
		if (improvement)
		{
		replace:
			Sol remember = sol;
			// execute replacement
			sol.replace_vertex(*besttour, candidate, position);
			if ((besttour->breakindex == int(besttour->seq.size()) - 1) && (ins->t[besttour->index].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "break pulled" << endl;
				pull_break(sol, besttour->index);
			}
			//sol.check();
			//cout << "hier" << endl;
			succes = true;
		}//end if improvement
	}//end while improvement
	return succes;
}

bool Moves::one_one_replace_nb(Sol& sol, int mode)
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestratio = 0.0;
		int best_pos_rem = -1;
		int best_pos_ins = -1;
		Ins::Vertex* bestcandidate = NULL;
		Sol::Tour* besttour = NULL;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			int t = sol.tourindex[d];//pick a random tour
			int endh = (int)sol.tours[sol.tourindex[d]].seq.size();
			for (int h = 1; h < endh - 1; ++h)// for all  inlcuded regular vertices in sol
			{
				//for all existing regular member vertices: remove 1 and revaluate solution
				Tour tourrem = sol.tours[sol.tourindex[d]];
				Ins::Vertex* r = tourrem.seq[h];//to be removed vertex
				tourrem.remove_vertex(h);
				int endj = (int)tourrem.seq.size();
				for (int j = 0; j < endj - 1; ++j)// for positions in tourrem
				{
					Ins::Vertex* x = tourrem.seq[j];//predecessor y
					Ins::Vertex* z = tourrem.seq[j + 1];//successor y
					int breakz = tourrem.action[j + 1];
					int nb_size = (int)x->nb[t].size();
					for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
					{
						Ins::Vertex* y = x->nb[t][i];//potential insertion at position j
						if ((sol.available[y->index] * y->score > r->score) && (x->nbi[t][y->index]) && (y->nbi[t][z->index]))//availability & score increase
						{
							if ((tourrem.weight + y->weight <= ins->t[t].W_max) && (tourrem.volume + y->volume <= ins->t[t].V_max))//cap constraint check
							{
								//gather departure time
								double currenttime = tourrem.deptime[j] + ins->t[t].EDT;//service bij x zit hier al in
								//travel time from x to y
								double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
								if (arrivaltime < y->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
								{
									arrivaltime = y->LTW[t];
								}
								if (arrivaltime > y->UTW[t])
								{
									continue;//infeasible
								}
								arrivaltime += y->serv;
								//travel time from y to z
								arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
								if (arrivaltime + breakz * (ins->breakdur) < z->LTW[t])
								{
									arrivaltime = z->LTW[t] - (breakz * ins->breakdur);
								}
								arrivaltime += z->serv + breakz * ins->breakdur;
								double shift = (arrivaltime - ins->t[t].EDT) - tourrem.deptime[j + 1];//increase in travel time
								if (shift <= tourrem.max_shift[j + 1])//check of het punt geinsert kan worden
								{
									improvement = true;
									double ratio = y->score - r->score;
									if (ratio > bestratio)
									{
										bestratio = ratio;
										besttour = &sol.tours[sol.tourindex[d]];
										best_pos_rem = h;
										best_pos_ins = j;//index after removal
										bestcandidate = y;
										if (mode == 0)//first improvement, otherwise best improvement
										{
											goto one_one_replace;
										}
									}
								}
							}//end cap constraints
						}//end availability &score check
					}//for all nb
				}//end for all positions in tourrem
			}//end for all included regular vertices in sol
		}//end for all tours
		if (improvement)
		{
			//delete first & insert after
			one_one_replace:
			sol.remove_vertex(*besttour, best_pos_rem);
			sol.insert_vertex(*besttour, bestcandidate, best_pos_ins);
			succes = true;
		}
	}//end while improvement
	return succes;
}//end one_one_replace_nb


bool Moves::two_one_replace_nb(Sol& sol, int mode)
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestratio = 0.0;
		int best_pos_rem1 = -1;
		int best_pos_rem2 = -1;
		int best_pos_ins = -1;
		Ins::Vertex* bestcandidate = NULL;
		Sol::Tour* besttour = NULL;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			int t = sol.tourindex[d];//pick a random tour
			int end = (int)sol.tours[sol.tourindex[d]].seq.size();
			for (int g = 1; g < end - 1; ++g)// for all  inlcuded regular vertices in sol
			{
				int end = (int)sol.tours[sol.tourindex[d]].seq.size();
				for (int h = 1; h < end - 1; ++h)// for all  inlcuded regular vertices in sol
				{
					if (g != h)
					{
						//for all existing regular member vertices: remove 1 and revaluate solution
						
						Sol::Tour tourrem = sol.tours[sol.tourindex[d]];
						Ins::Vertex* r = tourrem.seq[g];//to be removed vertex 1
						Ins::Vertex* s = tourrem.seq[h];//to be removed vertex 2
						int lostscore = r->score + s->score;
						tourrem.remove_vertex(g);
						if (g < h)
						{
							tourrem.remove_vertex(h - 1);
						}
						else
						{
							tourrem.remove_vertex(h);
						}
						int endj = (int)tourrem.seq.size();
						for (int j = 0; j < endj - 1; ++j)// for positions in tourrem
						{
							Ins::Vertex* x = tourrem.seq[j];//predecessor y
							Ins::Vertex* z = tourrem.seq[j + 1];//successor y
							int breakz = tourrem.action[j + 1];
							int nb_size = (int)x->nb[t].size();
							for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
							{
								Ins::Vertex* y = x->nb[t][i];//potential insertion at position j
								if ((sol.available[y->index] * y->score > lostscore) && (x->nbi[t][y->index]) && (y->nbi[t][z->index]))//availability & score increase
								{
									if ((tourrem.weight + y->weight <= ins->t[t].W_max) && (tourrem.volume + y->volume <= ins->t[t].V_max))//cap constraint check
									{
										//gather departure time
										double currenttime = tourrem.deptime[j] + ins->t[t].EDT;//service bij x zit hier al in
										//travel time from x to y
										double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
										if (arrivaltime < y->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
										{
											arrivaltime = y->LTW[t];
										}
										if (arrivaltime > y->UTW[t])
										{
											continue;//infeasible
										}
										arrivaltime += y->serv;
										//travel time from y to z
										arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
										if (arrivaltime + breakz * (ins->breakdur) < z->LTW[t])
										{
											arrivaltime = z->LTW[t] - (breakz * ins->breakdur);
										}
										arrivaltime += z->serv + breakz * ins->breakdur;
										double shift = (arrivaltime - ins->t[t].EDT) - tourrem.deptime[j + 1];//increase in travel time
										if (shift <= tourrem.max_shift[j + 1])//check of het punt geinsert kan worden
										{
											improvement = true;
											double ratio = y->score - lostscore;
											if (ratio > bestratio)
											{
												bestratio = ratio;
												besttour = &sol.tours[sol.tourindex[d]];
												best_pos_rem1 = g;
												best_pos_rem2 = h;
												best_pos_ins = j;//index after removal
												bestcandidate = y;
												if (mode == 0)//first improvement, otherwise best improvement
												{
													goto two_one_replace;
												}
											}
										}
									}//end cap constraints
								}//end availability &score check
							}//for all nb
						}//end for all positions in tourrem
					}//end if g & h are different
				}//end for all included regular vertices in sol
			}//end for all included regular vertices in sol
		}//end for all tours
		if (improvement)
		{
			two_one_replace:
			sol.remove_vertex(*besttour, best_pos_rem1);
			if (best_pos_rem1 < best_pos_rem2)
			{
				sol.remove_vertex(*besttour,best_pos_rem2-1);
			}
			else
			{
				sol.remove_vertex(*besttour,best_pos_rem2);
			}
			sol.insert_vertex(*besttour, bestcandidate, best_pos_ins);
			succes = true;
		}
	}//end while improvement
	return succes;
}//end two_one_replace_nb

bool Moves::one_two_replace(Sol& sol, int mode)
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestratio = 0.0;
		int best_pos_rem = -1;
		int best_pos_ins1 = -1;
		int best_pos_ins2 = -1;
		Ins::Vertex* bestcandidate1 = NULL;
		Ins::Vertex* bestcandidate2 = NULL;
		Sol::Tour* besttour = NULL;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			int t = sol.tourindex[d];//pick a random tour
			int endh = (int)sol.tours[sol.tourindex[d]].seq.size();
			for (int h = 1; h < endh - 1; ++h)// for all  inlcuded regular vertices in sol
			{
				//for all existing regular member vertices: remove 1 and revaluate solution
				Sol::Tour tourrem = sol.tours[sol.tourindex[d]];
				Ins::Vertex* r = tourrem.seq[h];//to be removed vertex
				tourrem.remove_vertex(h);
				//make candidate list
				vector<pair<Ins::Vertex*, Ins::Vertex*>> candidatelist;
				for (int i = 1; i < ins->maxvertices - 1; ++i)//for all regular vertices
				{
					for (int j = 1; j < ins->maxvertices - 1; ++j)//for all regular vertices
					{
						if (i != j)
						{
							Ins::Vertex* a = &ins->v[i];//candidate 1
							Ins::Vertex* b = &ins->v[j];//candidate 2
							if ((sol.available[a->index])&&(sol.available[b->index])&&(a->score+b->score>r->score))
							{//availability & score check
								if ((tourrem.weight + a->weight+b->weight <= ins->t[t].W_max) && (tourrem.volume + a->volume+b->volume <= ins->t[t].V_max))//cap constraint check
								{//capacity constraints
									candidatelist.push_back(pair<Ins::Vertex*, Ins::Vertex*>(a, b));
								}
							}
						}
					}
				}
				int candidatelistsize = (int) candidatelist.size();
				if (candidatelistsize > 0)
				{
					//check if both vertices can be inserted in the solution
					for (int c = 0; c < candidatelistsize - 1; ++c)//for all neighbours of the included vertex (non-enddepot)
					{
						int endj = (int)tourrem.seq.size();
						int pos1 = -1;
						int pos2 = -1;
						bool one = false;
						bool two = false;
						for (int j = 0; j < endj - 1; ++j)// for positions in tourrem
						{
							Ins::Vertex* x = tourrem.seq[j];//predecessor y
							Ins::Vertex* z = tourrem.seq[j + 1];//successor y
							int breakz = tourrem.action[j + 1];
							Ins::Vertex* y = candidatelist[c].first;//potential insertion at position j
							//gather departure time
							double currenttime = tourrem.deptime[j] + ins->t[t].EDT;//service bij x zit hier al in
							//travel time from x to y
							double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
							if (arrivaltime < y->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
							{
								arrivaltime = y->LTW[t];
							}
							if (arrivaltime > y->UTW[t])
							{
								continue;//infeasible
							}
							arrivaltime += y->serv;
							double arrivaltimey = arrivaltime;
							//travel time from y to z
							arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
							if (arrivaltime + breakz * (ins->breakdur) < z->LTW[t])
							{
								arrivaltime = z->LTW[t] - (breakz * ins->breakdur);
							}
							arrivaltime += z->serv + breakz * ins->breakdur;
							double shift = (arrivaltime - ins->t[t].EDT) - tourrem.deptime[j + 1];//increase in travel time
							if (shift <= tourrem.max_shift[j + 1])//check if first vertex can be inserted
							{
								one = true;
								pos1 = j;
								//check if second vertex can be inserted
								for (int k = 0; k < endj - 1; ++k)// for positions in tourrem
								{
									if (j != k)
									{
										Ins::Vertex* x = tourrem.seq[k];//predecessor y
										Ins::Vertex* z = tourrem.seq[k + 1];//successor y
										int breakz = tourrem.action[k + 1];
										Ins::Vertex* y = candidatelist[c].second;//potential insertion at position j
										//gather departure time
										double currenttime = tourrem.deptime[k] + ins->t[t].EDT;//service bij x zit hier al in
										//travel time from x to y
										double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
										if (arrivaltime < y->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
										{
											arrivaltime = y->LTW[t];
										}
										if (arrivaltime > y->UTW[t])
										{
											continue;//infeasible
										}
										arrivaltime += y->serv;
										//travel time from y to z
										arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
										if (arrivaltime + breakz * (ins->breakdur) < z->LTW[t])
										{
											arrivaltime = z->LTW[t] - (breakz * ins->breakdur);
										}
										arrivaltime += z->serv + breakz * ins->breakdur;
										double shift = (arrivaltime - ins->t[t].EDT) - tourrem.deptime[k + 1];//increase in travel time
										if (shift <= tourrem.max_shift[k + 1])//check of het punt geinsert kan worden
										{
											one = true;
											pos2 = k;
										}//end if shift
									}
									else
									{//special case were second vertex is inserted after first vertex on position j
										currenttime = arrivaltimey;
										x = y;
										y = candidatelist[c].second;
										double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
										if (arrivaltime < y->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
										{
											arrivaltime = y->LTW[t];
										}
										if (arrivaltime > y->UTW[t])
										{
											continue;//infeasible
										}
										arrivaltime += y->serv;
										//travel time from y to z
										arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
										if (arrivaltime + breakz * (ins->breakdur) < z->LTW[t])
										{
											arrivaltime = z->LTW[t] - (breakz * ins->breakdur);
										}
										arrivaltime += z->serv + breakz * ins->breakdur;
										double shift = (arrivaltime - ins->t[t].EDT) - tourrem.deptime[j + 1];//increase in travel time
										if (shift <= tourrem.max_shift[k + 1])//check of het punt geinsert kan worden
										{
											two = true;
											pos2 = k;
										}//end if shift
									}
									if ((one) && (two))
									{
										improvement = true;
										double ratio = (candidatelist[c].first->score + candidatelist[c].second->score) - r->score;
										if (ratio > bestratio)
										{
											bestratio = ratio;
											besttour = &sol.tours[sol.tourindex[d]];
											best_pos_rem = h;
											best_pos_ins1 = pos1;//index after removal
											best_pos_ins2 = pos2;//index after removal
											bestcandidate1 = candidatelist[c].first;
											bestcandidate2 = candidatelist[c].second;
											goto one_two_replace;
										}//end if
									}//end if both vertices can be inserted
								}//end for all positions in tourrem 2
							}//end if shift succes first insertion
						}//end for all positions in tourrem 1
					}//end for all pairs
				}//end if candidate list is not empty
			}//end for all included regular vertices in sol
		}//end for all tours
		if (improvement)
		{
			one_two_replace:
			sol.remove_vertex(*besttour, best_pos_rem);
			sol.insert_vertex(*besttour,bestcandidate1, best_pos_ins1);
			sol.insert_vertex(*besttour,bestcandidate2, best_pos_ins2+1);
			//sol.check();
			//cout << "debug here" << endl;
		}
	}//end while improvement
	return succes;
}//end one_two_replace_nb


boost::heap::priority_queue<One_one_rep_nb> Moves::one_one_replace_gen_nb(Sol& sol, int limit)
{
	boost::heap::priority_queue<One_one_rep_nb> adm_nb;
	int bestscore = -INT_MAX;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		//INSERT PART
		Sol::Tour& tour = sol.tours[d];
		int endj = (int)tour.seq.size();
		for (int j = 0; j < endj - 1; ++j)// for positions in the tour
		{
			Ins::Vertex* x = tour.seq[j];//predecessor y
			Ins::Vertex* z = tour.seq[j + 1];//successor y
			int breakz = tour.action[j + 1];//break on z
			int nb_size = (int)x->nb[d].size();
			for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
			{
				Ins::Vertex* y = x->nb[d][i];//potential insertion at position j
				if ((sol.available[y->index]) && (x->nbi[d][y->index]) && (y->nbi[d][z->index]))//availability & improvement to current best admissable
				{
					if ((y->score > bestscore) || (adm_nb.size() < limit))
					{
						if ((tour.weight + y->weight <= ins->t[d].W_max) && (tour.volume + y->volume <= ins->t[d].V_max))//cap constraint check
						{
							//gather departure time
							double currenttime = tour.deptime[j] + ins->t[d].EDT;//service bij x zit hier al in
							//travel time from x to y
							double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
							if (arrivaltime < y->LTW[d])//break inserten kan niet dus break kan ltw niet dichter brengen
							{
								arrivaltime = y->LTW[d];
							}
							if (arrivaltime > y->UTW[d])
							{
								continue;//infeasible
							}
							arrivaltime += y->serv;
							//travel time from y to z
							arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
							if (arrivaltime + breakz * (ins->breakdur) < z->LTW[d])
							{
								arrivaltime = z->LTW[d] - (breakz * ins->breakdur);
							}
							arrivaltime += z->serv + breakz * ins->breakdur;
							double shift = (arrivaltime - ins->t[d].EDT) - tour.deptime[j + 1];//increase in travel time
							if (shift <= tour.max_shift[j + 1])//check of het punt geinsert kan worden
							{
								adm_nb.push(One_one_rep_nb(d, -1, j, y, y->score));
								if (y->score > bestscore)
								{
									bestscore = y->score;
								}
							}
						}//end cap constraints
					}//end score or heap min limit
				}//end availability
			}//for all nb
		}//end for all positions in sol
		//REPLACE PART
		int endh = (int)sol.tours[d].seq.size();
		for (int h = 1; h < endh - 1; ++h)// for all  included regular vertices
		{
			//with replacing: for all existing regular member vertices: remove 1 and revaluate solution
			Tour tourrem = sol.tours[d];
			Ins::Vertex* r = tourrem.seq[h];//to be removed vertex
			tourrem.remove_vertex(h);
			int endj = (int)tourrem.seq.size();
			for (int j = 0; j < endj - 1; ++j)// for positions in tourrem
			{
				Ins::Vertex* x = tourrem.seq[j];//predecessor y
				Ins::Vertex* z = tourrem.seq[j + 1];//successor y
				int breakz = tourrem.action[j + 1];
				int nb_size = (int)x->nb[d].size();
				for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
				{
					Ins::Vertex* y = x->nb[d][i];//potential insertion at position j
					if ((sol.available[y->index]) && (x->nbi[d][y->index]) && (y->nbi[d][z->index]))//availability & improvement to current best admissable
					{
						if ((y->score-r->score > bestscore) || (adm_nb.size() < limit))
						{
							if ((tourrem.weight + y->weight <= ins->t[d].W_max) && (tourrem.volume + y->volume <= ins->t[d].V_max))//cap constraint check
							{
								//gather departure time
								double currenttime = tourrem.deptime[j] + ins->t[d].EDT;//service bij x zit hier al in
								//travel time from x to y
								double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
								if (arrivaltime < y->LTW[d])//break inserten kan niet dus break kan ltw niet dichter brengen
								{
									arrivaltime = y->LTW[d];
								}
								if (arrivaltime > y->UTW[d])
								{
									continue;//infeasible
								}
								arrivaltime += y->serv;
								//travel time from y to z
								arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
								if (arrivaltime + breakz * (ins->breakdur) < z->LTW[d])
								{
									arrivaltime = z->LTW[d] - (breakz * ins->breakdur);
								}
								arrivaltime += z->serv + breakz * ins->breakdur;
								double shift = (arrivaltime - ins->t[d].EDT) - tourrem.deptime[j + 1];//increase in travel time
								if (shift <= tourrem.max_shift[j + 1])//check of het punt geinsert kan worden
								{
									adm_nb.push(One_one_rep_nb(d, h, j, y, y->score - r->score));
									if (y->score - r->score > bestscore)
									{
										bestscore = y->score - r->score;
									}
								}
							}//end cap constraints
						}
					}//end availability 
				}//for all nb
			}//end for all positions in tourrem
		}//end for all included regular vertices in sol
	}//end for all tours
	return adm_nb;
}//end one_one_replace_gen_nb(Sol& sol)

boost::heap::priority_queue<One_two_rep_nb> Moves::one_two_replace_gen_nb(Sol& sol, int limit)
{
	int bestscore = -INT_MAX;
	boost::heap::priority_queue<One_two_rep_nb> adm_nb;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int endh = (int)sol.tours[d].seq.size();
		for (int h = 1; h < endh - 1; ++h)// for all included regular vertices in sol
		{
			//for all existing regular member vertices: remove 1 and revaluate solution
			
			Sol::Tour tourrem = sol.tours[d];
			Ins::Vertex* r = tourrem.seq[h];//to be removed vertex
			tourrem.remove_vertex(h);
			//make candidate list
			vector<pair<Ins::Vertex*, Ins::Vertex*>> candidatelist;
			for (int i = 1; i < ins->maxvertices - 1; ++i)//for all regular vertices
			{
				for (int j = 1; j < ins->maxvertices - 1; ++j)//for all regular vertices
				{
					if (i != j)
					{
						Ins::Vertex* a = &ins->v[i];//candidate 1
						Ins::Vertex* b = &ins->v[j];//candidate 2
						if ((sol.available[a->index]) && (sol.available[b->index]))
						{//availability
							
								if ((tourrem.weight + a->weight + b->weight <= ins->t[d].W_max) && (tourrem.volume + a->volume + b->volume <= ins->t[d].V_max))//cap constraint check
								{//capacity constraints
									candidatelist.push_back(pair<Ins::Vertex*, Ins::Vertex*>(a, b));
								}
						}
					}
				}
			}
			int candidatelistsize = (int)candidatelist.size();
			if (candidatelistsize > 0)
			{
				//check if both vertices can be inserted in the solution
				for (int c = 0; c < candidatelistsize - 1; ++c)//for all neighbours of the included vertex (non-enddepot)
				{
					if (((candidatelist[c].first->score + candidatelist[c].second->score) - r->score > bestscore) || (adm_nb.size() < limit))
					{
						int endj = (int)tourrem.seq.size();
						int pos1 = -1;
						int pos2 = -1;
						bool one = false;
						bool two = false;
						for (int j = 0; j < endj - 1; ++j)// for positions in tourrem
						{
							Ins::Vertex* x = tourrem.seq[j];//predecessor y
							Ins::Vertex* z = tourrem.seq[j + 1];//successor y
							int breakz = tourrem.action[j + 1];
							Ins::Vertex* y = candidatelist[c].first;//potential insertion at position j
							//gather departure time
							double currenttime = tourrem.deptime[j] + ins->t[d].EDT;//service bij x zit hier al in
							//travel time from x to y
							double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
							if (arrivaltime < y->LTW[d])//break inserten kan niet dus break kan ltw niet dichter brengen
							{
								arrivaltime = y->LTW[d];
							}
							if (arrivaltime > y->UTW[d])
							{
								continue;//infeasible
							}
							arrivaltime += y->serv;
							double arrivaltimey = arrivaltime;
							//travel time from y to z
							arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
							if (arrivaltime + breakz * (ins->breakdur) < z->LTW[d])
							{
								arrivaltime = z->LTW[d] - (breakz * ins->breakdur);
							}
							arrivaltime += z->serv + breakz * ins->breakdur;
							double shift = (arrivaltime - ins->t[d].EDT) - tourrem.deptime[j + 1];//increase in travel time
							if (shift <= tourrem.max_shift[j + 1])//check if first vertex can be inserted
							{
								one = true;
								pos1 = j;
								//check if second vertex can be inserted
								for (int k = 0; k < endj - 1; ++k)// for positions in tourrem
								{
									if (j != k)
									{
										Ins::Vertex* x = tourrem.seq[k];//predecessor y
										Ins::Vertex* z = tourrem.seq[k + 1];//successor y
										int breakz = tourrem.action[k + 1];
										Ins::Vertex* y = candidatelist[c].second;//potential insertion at position j
										//gather departure time
										double currenttime = tourrem.deptime[k] + ins->t[d].EDT;//service bij x zit hier al in
										//travel time from x to y
										double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
										if (arrivaltime < y->LTW[d])//break inserten kan niet dus break kan ltw niet dichter brengen
										{
											arrivaltime = y->LTW[d];
										}
										if (arrivaltime > y->UTW[d])
										{
											continue;//infeasible
										}
										arrivaltime += y->serv;
										//travel time from y to z
										arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
										if (arrivaltime + breakz * (ins->breakdur) < z->LTW[d])
										{
											arrivaltime = z->LTW[d] - (breakz * ins->breakdur);
										}
										arrivaltime += z->serv + breakz * ins->breakdur;
										double shift = (arrivaltime - ins->t[d].EDT) - tourrem.deptime[k + 1];//increase in travel time
										if (shift <= tourrem.max_shift[k + 1])//check of het punt geinsert kan worden
										{
											one = true;
											pos2 = k;
										}//end if shift
									}
									else
									{//special case were second vertex is inserted after first vertex on position j
										currenttime = arrivaltimey;
										x = y;
										y = candidatelist[c].second;
										double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
										if (arrivaltime < y->LTW[d])//break inserten kan niet dus break kan ltw niet dichter brengen
										{
											arrivaltime = y->LTW[d];
										}
										if (arrivaltime > y->UTW[d])
										{
											continue;//infeasible
										}
										arrivaltime += y->serv;
										//travel time from y to z
										arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
										if (arrivaltime + breakz * (ins->breakdur) < z->LTW[d])
										{
											arrivaltime = z->LTW[d] - (breakz * ins->breakdur);
										}
										arrivaltime += z->serv + breakz * ins->breakdur;
										double shift = (arrivaltime - ins->t[d].EDT) - tourrem.deptime[j + 1];//increase in travel time
										if (shift <= tourrem.max_shift[k + 1])//check of het punt geinsert kan worden
										{
											two = true;
											pos2 = k;
										}//end if shift
									}
									if ((one) && (two))
									{
										adm_nb.push(One_two_rep_nb(d, h, pos1, pos2, candidatelist[c].first, candidatelist[c].second, (candidatelist[c].first->score + candidatelist[c].second->score) - r->score));
										if ((candidatelist[c].first->score + candidatelist[c].second->score) - r->score > bestscore)
										{
											bestscore = (candidatelist[c].first->score + candidatelist[c].second->score) - r->score;
										}
									}//end if both vertices can be inserted
								}//end for all positions in tourrem 2
							}//end if shift succes first insertion
						}//end for all positions in tourrem 1
					}//end score & nb limit check
				}//end for all pairs
			}//end if candidate list is not empty
		}//end for all included regular vertices in sol
	}//end for all tours
	return adm_nb;
}// end one_two_replace_gen_nb(Sol& sol)


boost::heap::priority_queue<Two_one_rep_nb> Moves::two_one_replace_gen_nb(Sol& sol, int limit)
{
	int bestscore = -INT_MAX;
	boost::heap::priority_queue<Two_one_rep_nb> adm_nb;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int end = (int)sol.tours[d].seq.size();
		for (int g = 1; g < end - 1; ++g)// for all  included regular vertices in sol
		{
			for (int h = 1; h < end - 1; ++h)// for all  inlcuded regular vertices in sol
			{
				if (g != h)
				{
					//for all existing regular member vertices: remove 1 and revaluate solution
					Sol::Tour tourrem = sol.tours[d];
					Ins::Vertex* r = tourrem.seq[g];//to be removed vertex 1
					Ins::Vertex* s = tourrem.seq[h];//to be removed vertex 2
					int lostscore = r->score + s->score;
					tourrem.remove_vertices(g, h);
					int endj = (int)tourrem.seq.size();
					for (int j = 0; j < endj - 1; ++j)// for positions in tourrem
					{
						Ins::Vertex* x = tourrem.seq[j];//predecessor y
						Ins::Vertex* z = tourrem.seq[j + 1];//successor y
						int breakz = tourrem.action[j + 1];
						int nb_size = (int)x->nb[d].size();
						for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
						{
							Ins::Vertex* y = x->nb[d][i];//potential insertion at position j
							if ((sol.available[y->index]) && (x->nbi[d][y->index]) && (y->nbi[d][z->index]))//availability
							{
								if ((y->score - lostscore > bestscore) || (adm_nb.size() < limit))
								{
									if ((tourrem.weight + y->weight <= ins->t[d].W_max) && (tourrem.volume + y->volume <= ins->t[d].V_max))//cap constraint check
									{
										//gather departure time
										double currenttime = tourrem.deptime[j] + ins->t[d].EDT;//service bij x zit hier al in
										//travel time from x to y
										double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
										if (arrivaltime < y->LTW[d])//break inserten kan niet dus break kan ltw niet dichter brengen
										{
											arrivaltime = y->LTW[d];
										}
										if (arrivaltime > y->UTW[d])
										{
											continue;//infeasible
										}
										arrivaltime += y->serv;
										//travel time from y to z
										arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
										if (arrivaltime + breakz * (ins->breakdur) < z->LTW[d])
										{
											arrivaltime = z->LTW[d] - (breakz * ins->breakdur);
										}
										arrivaltime += z->serv + breakz * ins->breakdur;
										double shift = (arrivaltime - ins->t[d].EDT) - tourrem.deptime[j + 1];//increase in travel time
										if (shift <= tourrem.max_shift[j + 1])//check of het punt geinsert kan worden
										{
											adm_nb.push(Two_one_rep_nb(d, g, h, j, y, y->score - lostscore));
											if (y->score - lostscore > bestscore)
											{
												bestscore = y->score - lostscore;
											}
										}
									}//end cap constraints
								}//end score & nb limit check
							}//end availability/nb check
						}//for all nb
					}//end for all positions in tourrem
				}//end if g & h are different
			}//end for all included regular vertices in sol
		}//end for all included regular vertices in sol
	}//end for all tours
	return adm_nb;
}// end two_one_replace_gen_nb(Sol& sol)


void Moves::pull_break(Sol& sol, int t)
{
	int end = (int)sol.tours[t].seq.size();
	Sol::Tour* tour = &sol.tours[t];
	int earliestbreakindex = end - 1;
	for (int i = 0; i < end; ++i)
	{
		Ins::Vertex* p = tour->seq[i];
		int breakcurrent = tour->action[i];
		double arrivaltime = (ins->t[t].EDT + tour->deptime[i]) - (p->serv + breakcurrent * ins->breakdur);
		if(arrivaltime >= ins->breakstart)//check if arrivaltime
		{
			earliestbreakindex = i;
			break;
		}
	}
	for (int i = earliestbreakindex; i < tour->breakindex; ++i)//huidige breakindex niet evalueren die ken je al
	{
		double currenttime = ins->t[t].EDT + tour->deptime[i - 1];
		bool feasible = true;
		for (int j = i - 1; j < end - 1; ++j)
		{
			Ins::Vertex* last = tour->seq[j];
			Ins::Vertex* current = tour->seq[j + 1];
			int breakcurrent = 0;
			if (j+1 == i)
			{
				breakcurrent = 1;
			}
			double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
			//you can arrive too early to take a break
			if (arrivaltime < ins->breakstart)
			{
				feasible = false;
				break;
			}
			arrivaltime += breakcurrent * ins->breakdur;
			if (arrivaltime < current->LTW[t])
			{
				arrivaltime = current->LTW[t];
			}
			//utw check want je kan later aankomen door de break vroeger te schedulen
			if (arrivaltime > current->UTW[t])
			{
				feasible = false;
				break;
			}
			arrivaltime += current->serv;
			currenttime = arrivaltime;
		}
		if ((feasible) && (currenttime <= ins->t[t].EDT + tour->deptime.back()))
		{
			tour->action[tour->breakindex] = 0;
			tour->action[i] = 1;
			tour->breakindex = i;
			currenttime = ins->t[t].EDT + tour->deptime[i - 1];
			for (int j = i - 1; j < end - 1; ++j)
			{
				Ins::Vertex* last = tour->seq[j];
				Ins::Vertex* current = tour->seq[j + 1];
				int breakcurrent = tour->action[j + 1];
				double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
				arrivaltime += breakcurrent * ins->breakdur;
				if (arrivaltime < current->LTW[t])
				{
					arrivaltime = current->LTW[t];
				}
				arrivaltime += current->serv;
				tour->deptime[j + 1] = arrivaltime - ins->t[t].EDT;
				currenttime = arrivaltime;
			}
			tour->calc_maxshift();
			break;
		}//end for improvement
	}//end for all break combinations
}

void Moves::reschedule_breaks(Sol& sol)
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		Sol::Tour* tour = &sol.tours[d];
		int end = (int) tour->seq.size();
		int earliestbreakindex = end - 1;
		int bestindex = -1;
		double bestdecrease = 0;
		bool improvement = false;
		for (int i = 0; i < end; ++i)
		{
			Ins::Vertex* p = tour->seq[i];
			int breakcurrent = tour->action[i];
			if ((ins->t[d].EDT + tour->deptime[i]) - (p->serv+breakcurrent*ins->breakdur) >= ins->breakstart)
			{
				earliestbreakindex = i;
				break;
			}
		}
		for (int i = earliestbreakindex; i < tour->breakindex; ++i)//huidige breakindex niet evalueren die ken je al
		{
			double currenttime = ins->t[d].EDT + tour->deptime[i - 1];
			bool feasible = true;
			for (int j = i - 1; j < end - 1; ++j)
			{
				Ins::Vertex* last = tour->seq[j];
				Ins::Vertex* current = tour->seq[j + 1];
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
			if ((feasible) && (currenttime < ins->t[d].EDT + tour->deptime.back()))
			{
				double decrease = (ins->t[d].EDT + tour->deptime.back()) - currenttime;
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
			tour->action[tour->breakindex] = 0;
			tour->action[bestindex] = 1;
			tour->breakindex = bestindex;
			double currenttime = ins->t[d].EDT + tour->deptime[bestindex - 1];
			for (int j = bestindex - 1; j < end - 1; ++j)
			{
				Ins::Vertex* last = tour->seq[j];
				Ins::Vertex* current = tour->seq[j + 1];
				int breakcurrent = tour->action[j + 1];
				double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
				if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[d])
				{
					arrivaltime = current->LTW[d] - (breakcurrent * ins->breakdur);
				}
				arrivaltime += current->serv+breakcurrent*ins->breakdur;
				tour->deptime[j + 1] = arrivaltime - ins->t[d].EDT;
				currenttime = arrivaltime;
			}
			tour->max_shift.back() = (ins->t[d].T_max - tour->deptime.back());
			double departuretime = 0;
			double arrivaltime = ins->t[d].LAT;
			if (tour->action.back() == 1)//break op enddepot
			{
				if (ins->t[d].LAT > ins->breakend + ins->breakdur)
				{
					//cout<<"path: "<<d<< "pull break break op enddepot verhindert een maxshift: " << endl;
					tour->max_shift.back() = (ins->breakend + ins->breakdur) - (tour->deptime.back() + ins->t[d].EDT);
					arrivaltime = ins->breakend;//zoals hieronder service of enkel break in dit geval ervan aftrekken
				}
			}
			Ins::Vertex* y;
			Ins::Vertex* z;
			for (int i = 0; i < end - 2; ++i)// depots don't count
			{
				//define the 2 elements
				y = tour->seq[end - (i + 2)];
				int breaki = tour->action[end - (i + 2)];
				z = tour->seq[end - (i + 1)];
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
				tour->max_shift[end - (i + 2)] = departuretime - (tour->deptime[end - (i + 2)] + ins->t[d].EDT);
				//reset variable for the calculation of next point
				arrivaltime = departuretime - y->serv+breaki*ins->breakdur;
			}// end for i
		}//end if improvement
	}//end for all tours
}//end reschedule_breaks

bool Moves::exchange2_nb(Sol& sol)//exchange two vertices between two routes
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 0;
		Sol::Tour*  bestpatha=NULL;
		Sol::Tour* bestpathb=NULL;
		int bestindexa;
		int bestindexb;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			Sol::Tour* tourd = &sol.tours[d];
			for (int i = 0; i < (int)tourd->seq.size() - 2; ++i)
			{
				Ins::Vertex* x = tourd->seq[i + 1];
				Ins::Vertex* y = tourd->seq[i + 2];
				double oldtraveltime = (tourd->deptime[i + 2] - (y->serv + tourd->action[i + 2] * ins->breakdur + x->serv + tourd->action[i + 1] * ins->breakdur)) - tourd->deptime[i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					Sol::Tour* toure= &sol.tours[e];
					if (d != e)
					{
						for (int j = 0; j < (int) toure->seq.size() - 1; ++j)
						{
							Ins::Vertex* a = toure->seq[j];
							Ins::Vertex* b = toure->seq[j + 1];
							double departuretime = ins->t[e].EDT + toure->deptime[j];
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
								arrivaltime += b->serv + toure->action[j + 1] * ins->breakdur;
								double shiftmaster = (arrivaltime - ins->t[e].EDT) - toure->deptime[j + 1];//increase in travel time when inserting
								if ((shiftmaster <= toure->max_shift[j + 1]))//point can be inserted
								{
									improvement = true;
									bestpatha = toure;
									bestpathb = tourd;//slave
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
			Ins::Vertex* x = bestpathb->seq[bestindexb - 1];//point before insertion
			Ins::Vertex* exch = bestpathb->seq[bestindexb];//point that will be exchanged
			bestpathb->score -= exch->score;
			bestpathb->volume -= exch->volume;
			bestpathb->weight -= exch->weight;
			Ins::Vertex* candidate = NULL;
			int nb_size = (int)x->nb[bestpathb->index].size();
			for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex
			{
				Ins::Vertex* y = x->nb[bestpathb->index][i];//point that might be inserted
				Ins::Vertex* z = bestpathb->seq[bestindexb + 1];//point to shift
				int breakz = bestpathb->action[bestindexb + 1];
				if ((sol.available[y->index]) && (y->nbi[bestpathb->index][z->index]))
				{
					if ((bestpathb->volume + y->nb[bestpathb->index][i]->volume <= ins->t[bestpathb->index].V_max) && (bestpathb->weight + y->nb[bestpathb->index][i]->weight <= ins->t[bestpathb->index].W_max))
					{
						//gather departure time
						double currenttime = bestpathb->deptime[bestindexb - 1] + ins->t[bestpathb->index].EDT;//service bij x zit hier al in
						//travel time from x to y
						double arrivaltime = ins->arrival_time(x->con[y->index], currenttime);
						if (arrivaltime < y->LTW[bestpathb->index])
						{
							arrivaltime = y->LTW[bestpathb->index];
						}
						if (arrivaltime > y->UTW[bestpathb->index])
						{
							continue;//mag je al stoppen met rekenen voor dit punt op deze positie
						}
						arrivaltime += y->serv + bestpathb->action[bestindexb] * ins->breakdur;
						//travel time from y to z
						arrivaltime = ins->arrival_time(y->con[z->index], arrivaltime);
						if (arrivaltime < z->LTW[bestpathb->index])
						{
							arrivaltime = z->LTW[bestpathb->index];
						}
						arrivaltime += z->serv + breakz * ins->breakdur;
						double shift = (arrivaltime - ins->t[bestpathb->index].EDT) - bestpathb->deptime[bestindexb + 1];//increase in travel time
						if (shift <= bestpathb->max_shift[bestindexb + 1])//check of het punt geinsert kan worden
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
				bestpathb->seq[bestindexb] = candidate;
				bestpathb->score += candidate->score;// update score of the new solution
				sol.score += candidate->score;// update score of the new solution
				bestpathb->volume += candidate->volume;
				bestpathb->weight += candidate->weight;
				//update travel time solution vector 
				double currenttime = bestpathb->deptime[bestindexb - 1] + ins->t[bestpathb->index].EDT;
				for (int u = position; u < bestpathb->seq.size(); ++u)
				{
					//gather departure time and corresponding time slot
					Ins::Vertex* o = bestpathb->seq[u - 1];
					Ins::Vertex* p = bestpathb->seq[u];
					//travel time from van o to p
					double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
					if (arrivaltime < p->LTW[bestpathb->index])
					{
						arrivaltime = p->LTW[bestpathb->index];
					}
					arrivaltime += p->serv + (bestpathb->action[u] * ins->breakdur);
					bestpathb->max_shift[u] = (bestpathb->deptime[u] + bestpathb->max_shift[u]) - (arrivaltime - ins->t[bestpathb->index].EDT);
					//old arrival time + old maxshift - new arrival time = new max_shift
					bestpathb->deptime[u] = arrivaltime - ins->t[bestpathb->index].EDT;
					currenttime = arrivaltime;
				}//end for
				//the value of max_shift j+1 is now wrong but will be soon be updated
				//complete update of max_shift for the included vertices before the insertion
				double departuretime = 0;
				double arrivaltime = (bestpathb->deptime[position + 1] + ins->t[bestpathb->index].EDT + bestpathb->max_shift[position + 1]) - bestpathb->seq[position + 1]->serv + bestpathb->action[position + 1] * ins->breakdur;//service time eraftrekken
				for (int vz = position; vz > 0; --vz)
				{
					//define the 2 elements
					Ins::Vertex* f = bestpathb->seq[vz];
					int breakf = bestpathb->action[vz];
					Ins::Vertex* g = bestpathb->seq[vz + 1];
					//find proper departure time and time slot
					departuretime = ins->departure_time(f->con[g->index], arrivaltime);	//utw check
					if (departuretime > f->UTW[bestpathb->index] + f->serv + breakf * ins->breakdur)
					{
						departuretime = f->UTW[bestpathb->index] + f->serv + breakf * ins->breakdur;
					}
					if (breakf == 1)
					{
						if (departuretime > ins->breakend + ins->breakdur)
						{
							departuretime = ins->breakend + ins->breakdur;
						}
					}
					bestpathb->max_shift[vz] = departuretime - (bestpathb->deptime[vz] + ins->t[bestpathb->index].EDT);
					departuretime -= f->serv + breakf * ins->breakdur;
					arrivaltime = departuretime;
				}// end for vz

				//update master
				bestpatha->seq.insert(bestpatha->seq.begin() + bestindexa + 1, exch);//insert exch after x
				bestpatha->deptime.insert(bestpatha->deptime.begin() + bestindexa + 1, 0);//insert temporary value
				bestpatha->action.insert(bestpatha->action.begin() + bestindexa + 1, 0);//insert regular visit action change later when necessary
				if (bestindexa < bestpatha->breakindex)
				{
					bestpatha->breakindex += 1;//update breakindex
				}
				bestpatha->max_shift.insert(bestpatha->max_shift.begin() + bestindexa + 1, 0);
			    bestpatha->score += exch->score;// update score of the new solution
				//score of solution remains the same
				bestpatha->volume += exch->volume;
				bestpatha->weight += exch->weight;
				//update travel time solution vector 
				currenttime = bestpatha->deptime[bestindexa] + ins->t[bestpatha->index].EDT;
				for (int u = bestindexa + 1; u < bestpatha->seq.size(); ++u)
				{
					//gather departure time and corresponding time slot
					Ins::Vertex* o = bestpatha->seq[u - 1];
					Ins::Vertex* p =bestpatha->seq[u];
					//travel time from van o to p
					double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
					if (arrivaltime < p->LTW[bestpatha->index])
					{
						arrivaltime = p->LTW[bestpatha->index];
					}
					arrivaltime += p->serv + bestpatha->action[u] * ins->breakdur;
					bestpatha->max_shift[u] = (bestpatha->deptime[u] + bestpatha->max_shift[u]) - (arrivaltime - ins->t[bestpatha->index].EDT);
					//old arrival time + old maxshift - new arrival time = new max_shift
					bestpatha->deptime[u] = arrivaltime - ins->t[bestpatha->index].EDT;
					currenttime = arrivaltime;
				}//end for
				//the value of max_shift j+1 is now wrong but will be soon be updated
				//complete update of max_shift for the included vertices before the insertion
				departuretime = 0;
				arrivaltime = (bestpatha->deptime[bestindexa + 2] + ins->t[bestpatha->index].EDT + bestpatha->max_shift[bestindexa + 2]) - bestpatha->seq[bestindexa + 2]->serv + bestpatha->action[bestindexa + 2] * ins->breakdur;//service time eraftrekken
				for (int vz = bestindexa + 1; vz > 0; --vz)
				{
					//define the 2 elements
					Ins::Vertex* f = bestpatha->seq[vz];
					int breakf = bestpatha->action[vz];
					Ins::Vertex* g = bestpatha->seq[vz + 1];
					//find proper departure time and time slot
					departuretime = ins->departure_time(f->con[g->index], arrivaltime);	//utw check
					if (departuretime > f->UTW[bestpatha->index] + f->serv + breakf * ins->breakdur)
					{
						departuretime = f->UTW[bestpatha->index] + f->serv + breakf * ins->breakdur;
					}
					if (breakf == 1)
					{
						if (departuretime > ins->breakend + ins->breakdur)
						{
							departuretime = ins->breakend + ins->breakdur;
						}
					}
					bestpatha->max_shift[vz] = departuretime - (bestpatha->deptime[vz] + ins->t[bestpatha->index].EDT);
					departuretime -= f->serv + breakf * ins->breakdur;
					arrivaltime = departuretime;
				}// end for vz
				succes = true;
			}
			else
			{
				improvement = false;
				bestpathb->score += exch->score;
				bestpathb->volume += exch->volume;
				bestpathb->weight += exch->weight;
			}
			//sol.check();
		}//end if improvement
	}//end while improvement
	return succes;
}//end exchange


bool Moves::swap_nb(Sol& sol, int mode)
{
	bool succes = false;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		Sol::Tour& tour = sol.tours[d];
		bool improvement = true;
		while (improvement)
		{
			double bestdelta = 0.0001;
			int besti=-1;
			int bestj=-1;
			improvement = false;
			int end = (int)tour.seq.size();
			for (int i = 1; i < end - 2; ++i)// every vertex except depots and second last real vertex
			{
				for (int j = i + 1; j < end - 1; ++j)// every vertex except depots and second last real vertex
				{
					Ins::Vertex* k = tour.seq[i - 1];//pred swap partner 1
					Ins::Vertex* y = tour.seq[i];//swap partner 1
					Ins::Vertex* l = tour.seq[i + 1];//suc partner 1
					Ins::Vertex* m = tour.seq[j - 1];//pred swap partner 2
					Ins::Vertex* z = tour.seq[j];//swap partner 2
					Ins::Vertex* n = tour.seq[j + 1];//suc partner 2
					if (z == l)//fix pointers in special case y and z are directly linked
					{
						l = tour.seq[j - 1];//suc partner 1
						m = tour.seq[i + 1];//pred swap partner 2
					}
					if ((k->nbi[d][z->index]) && (m->nbi[d][y->index]) && (z->nbi[d][l->index]) && (y->nbi[d][n->index]))
					{//capacity constraints remain equal
						double delta_tt = 0;
						double departuretime = ins->t[d].EDT + tour.deptime[i - 1];//bij punt voor y
						double currenttime = departuretime;
						double newtraveltime;
						double oldtraveltime;
						//execute temporary swap using a auxilary linker variable break stays at the same position
						vector<int>linker(end, 0);
						for (int z = 0; z < end; ++z)
						{
							linker[z] = z;
						}
						linker[i] = j;
						linker[j] = i;
						for (int l = i - 1; l < j + 1; ++l)
						{
							Ins::Vertex* first = tour.seq[linker[l]];
							Ins::Vertex* second = tour.seq[linker[l + 1]];
							double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
							//do not account for break repositioning so don't use linker on action
							if (arrivaltime + tour.action[l + 1] * (ins->breakdur) < second->LTW[d])
							{
								arrivaltime = second->LTW[d] - (tour.action[l + 1] * ins->breakdur);
							}
							if (arrivaltime > second->UTW[d])
							{
								currenttime = DBL_MAX;
								break;
							}
							arrivaltime += second->serv + (tour.action[l + 1] * ins->breakdur);
							currenttime = arrivaltime;
						}// end l
						newtraveltime = currenttime - departuretime;
						oldtraveltime = (ins->t[d].EDT+tour.deptime[j+1])-departuretime;
						delta_tt = oldtraveltime-newtraveltime;
						if (delta_tt > bestdelta)//local evaluation
						{
							bool reqbreak = true;
							currenttime = ins->t[d].EDT;
							for (int l = 0; l < end - 1; ++l)//global evaluation
							{
								Ins::Vertex* first = tour.seq[linker[l]];
								Ins::Vertex* second = tour.seq[linker[l + 1]];//can be the end depot
								double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
								//account for break
								if ((reqbreak) && ((max(second->LTW[d] - ins->breakdur,arrivaltime) >= ins->breakstart) || (second->index == ins->maxvertices - 1)))
								{
									arrivaltime += ins->breakdur;
									reqbreak = false;
								}
								if (arrivaltime < second->LTW[d])
								{
									arrivaltime = second->LTW[d];
								}
								if (arrivaltime > second->UTW[d])
								{
									currenttime = DBL_MAX;
									break;
								}
								arrivaltime += second->serv;
								currenttime = arrivaltime;
							}// end l
							if (reqbreak == true)//als je aankomt bij het einddepot en nog steeds geen break genomen hebt
							{
								currenttime += ins->breakdur;//breaktime bijtellen bij aankomst tijd bij einddepot
							}
							newtraveltime = currenttime - departuretime;
							oldtraveltime = ((ins->t[d].EDT + tour.deptime[end - 1]) - departuretime);
							delta_tt = oldtraveltime-newtraveltime;
							if (delta_tt > bestdelta)//global evaluation
							{
								improvement = true;
								bestdelta = delta_tt;
								besti = i;
								bestj = j;
								if (mode == 0)
								{
									goto swap;
								}
							}//end global evaluation
						}// end local evaluation
					}//end neighbourhood check
				}//end j
			}// end i
			if (improvement)
			{
				swap:
				Sol remember = sol;
				tour.swap_vertices(besti, bestj);
				double actualdecrease = remember.tours[d].deptime.back() - tour.deptime.back();
				if (abs(bestdelta - actualdecrease) > 0.01)
				{
					cout << "error swap" << endl;
				}
				succes = true;
			}
		}//end while improvement
	}//end for all d
	return succes;
}//end swap

bool Moves::two_opt_nb(Sol& sol,int mode)
{
	bool succes = false;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		Sol::Tour& tour = sol.tours[d];
		bool improvement = true;
		while (improvement)
		{
			double bestdelta = 0.0001;
			int besti = -1;
			int bestj = -1;
			improvement = false;
			int end = (int)tour.seq.size();
			for (int i = 1; i < end - 1; ++i)// depots don't count
			{
				for (int j = i + 1; j < end - 1; ++j)
				{
					Ins::Vertex* k = tour.seq[i - 1];//pred partner 1
					Ins::Vertex* y = tour.seq[i];//partner 1
					Ins::Vertex* l = tour.seq[i + 1];//suc partner 1
					Ins::Vertex* m = tour.seq[j - 1];//pred  partner 2
					Ins::Vertex* z = tour.seq[j];// partner 2
					Ins::Vertex* n = tour.seq[j + 1];//suc partner 2
					if ((k->nbi[d][z->index]) && (z->nbi[d][m->index]) && (l->nbi[d][y->index]) && (y->nbi[d][n->index]))
					{
						//temporary opt and leave break on seq position
						vector<int>linker(end, 0);
						for (int z = 0; z < end; ++z)
						{
							linker[z] = z;
						}
						for (int f = 0; f < 1 + (j - i) / 2; ++f)
						{
							linker[j - f] = i + f;
							linker[i + f] = j - f;
						}
						
						bool reqbreak = false;
						if (i <= tour.breakindex)
						{//break is positioned after i, so might need to be repositioned
							reqbreak = true;
						}
						double delta_tt = 0;
						double departuretime = ins->t[d].EDT + tour.deptime[i - 1];
						double currenttime = departuretime;
						double arrivaltime;
						for (int l = i - 1; l < j + 1; ++l)
						{
							Ins::Vertex* first = tour.seq[linker[l]];
							Ins::Vertex* second = tour.seq[linker[l + 1]];
							arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
							if ((reqbreak) && ((max(second->LTW[d] - ins->breakdur,arrivaltime) >= ins->breakstart) || (second->index == ins->maxvertices - 1)))
							{
								arrivaltime += ins->breakdur;
								reqbreak = false;
							}
							if (arrivaltime < second->LTW[d])
							{
								arrivaltime = second->LTW[d];
							}
							if (arrivaltime > second->UTW[d])
							{
								currenttime = DBL_MAX;
								break;
							}
							currenttime = arrivaltime + second->serv;
						}// end l
						//int remember = sol.traveltime[d].back();
						delta_tt = (ins->t[d].EDT + tour.deptime[j + 1])-currenttime;
						if (delta_tt > 0)//local evaluation
						{
							Tour temptour = tour;
							temptour.opt_vertices(i,j);
							delta_tt = tour.deptime.back() - temptour.deptime.back();
							if (delta_tt>bestdelta)//global evaluation
							{
								improvement = true;
								bestdelta = delta_tt;
								besti = i;
								bestj = j;
							}
						}//end executed 2 opt
					}//end neighborhood check
				}// end for al j
			}// end for all i
			if (improvement)
			{
				Sol remember = sol;
				tour.opt_vertices(besti,bestj);
				double actualdecrease = remember.tours[d].deptime.back() - tour.deptime.back();
				if (abs(bestdelta - actualdecrease) > 0.01)
				{
					cout << "error 2opt" << endl;
				}
				succes = true;
			}
		}//end while improvement
	}//end for all paths
	return succes;
}

bool Moves::move_nb(Sol& sol,int mode)//move vertex x from tour d to tour e in order to save travel time
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 0.0001;
		Sol::Tour *bestd=NULL;
		Sol::Tour *beste=NULL;
		int besti;
		int bestj;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			Sol::Tour* tourd = &sol.tours[d];
			for (int i = 0; i < int(tourd->seq.size()) - 2; ++i)
			{
				Ins::Vertex* w = tourd->seq[i];
				Ins::Vertex* x = tourd->seq[i + 1];
				Ins::Vertex* y = tourd->seq[i + 2];
				int breaky = 0;
				double ttwxy = (tourd->deptime[i + 2] - (y->serv + tourd->action[i + 2] * ins->breakdur + x->serv + tourd->action[i + 1] * ins->breakdur)) - tourd->deptime[i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					Sol::Tour* toure = &sol.tours[e];
					if (d != e)//no relocate on the same path
					{
						for (int j = 0; j < int(toure->seq.size()) - 2; ++j)
						{
							Ins::Vertex* a = toure->seq[j];
							Ins::Vertex* b = toure->seq[j + 1];
							if ((a->nbi[e][x->index]) && (x->nbi[e][b->index]) && (w->nbi[d][y->index]))
							{
								if ((toure->weight + x->weight <= ins->t[e].W_max) && (toure->volume + x->volume <= ins->t[e].V_max))
								{
									//local evaluation on path e
									double ttab = (toure->deptime[j + 1] - (b->serv + toure->action[j + 1] * ins->breakdur)) - toure->deptime[j];
									//calculate axb
									double departuretime = ins->t[e].EDT + toure->deptime[j];
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
									if (arrivaltime + toure->action[j + 1] * ins->breakdur < b->LTW[e])
									{
										arrivaltime = b->LTW[e] - (toure->action[j + 1] * ins->breakdur);
									}
									if (arrivaltime > b->UTW[e])
									{
										continue;
									}
									ttaxb += arrivaltime - departuretime;
									double arrivalb = arrivaltime + b->serv + toure->action[j + 1] * ins->breakdur;

									//local evaluation on path d
									//due to potential traveltime decrease break can come to early
									bool reqbreak = false;
									if (tourd->breakindex >= i + 1)
									{
										reqbreak = true;
									}
									//calculate wy
									departuretime = ins->t[d].EDT + tourd->deptime[i];
									//traveltime w to y
									arrivaltime = ins->arrival_time(w->con[y->index], departuretime);
									//break niet meetellen voor local evaluation
									if ((reqbreak) && ((max(y->LTW[d]-ins->breakdur,arrivaltime) >= ins->breakstart) || (y->index == ins->maxvertices - 1)))
									{//break na y wordt bij global evaluation in rekening gebracht
										if (arrivaltime > ins->breakend)//pushing to break to right after removal does not work
										{
											continue;
										}
										reqbreak = false;
										breaky = 1;
										arrivaltime += ins->breakdur;
									}
									if (arrivaltime < y->LTW[d])
									{
										arrivaltime = y->LTW[d];
									}
									if (arrivaltime > y->UTW[d])
									{
										continue;
									}
									double ttwy = arrivaltime - departuretime;
									double arrivaly = arrivaltime + y->serv;
									double localincrease = (ttaxb - ttab) + x->serv;
									double localdecreasetotal = ttwxy + ttab - (ttwy + ttaxb);
									if ((localincrease <= toure->max_shift[j + 1]) && (localdecreasetotal > bestdecrease))//local improvement check
									{
										//check enddepot time on path d
										Tour temptourd = *tourd;
										temptourd.remove_vertex(i + 1);
										double globaldecreasetotal = (tourd->deptime.back() - temptourd.deptime.back());
										
										//check enddepot time on path e
										double currenttime = arrivalb;
										for (int m = j + 2; m < (int)toure->seq.size(); ++m)
										{
											//gather departure time and corresponding time slot
											Ins::Vertex* o = toure->seq[m - 1];
											Ins::Vertex* p = toure->seq[m];
											//travel time from van o to p
											double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
											if (arrivaltime + toure->action[m] * ins->breakdur < p->LTW[e])
											{
												arrivaltime = p->LTW[e] - (toure->action[m] * ins->breakdur);
											}
											arrivaltime += p->serv + toure->action[m] * ins->breakdur;
											currenttime = arrivaltime;
										}
										//cout << "enddepot time e: " << currenttime - t[e].EDT << endl;
										//cout << "increase e: " << (currenttime - t[e].EDT) - sol.traveltime[e].back() << endl;
										globaldecreasetotal -= ((currenttime - ins->t[e].EDT) - toure->deptime.back());
										if (globaldecreasetotal > bestdecrease)
										{
											bestdecrease = globaldecreasetotal;
											beste = toure;
											bestj = j;
											bestd = tourd;
											besti = i + 1;
											improvement = true;
											if (mode == 0)
											{
												goto move;
											}
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
			move:
			Sol remember = sol;
			Ins::Vertex* candidate = bestd->seq[besti];
			sol.remove_vertex(*bestd, besti);//remove vertex from path d
			sol.insert_vertex(*beste, candidate, bestj);//insert vertex on path e
			double actualdecrease = 0.0;
			for (int t = 0; t < (int)sol.tours.size(); ++t)
			{
				actualdecrease += remember.tours[t].deptime.back() - sol.tours[t].deptime.back();
			}
			if (abs(bestdecrease - actualdecrease) > 0.01)
			{
				cout << "error move" << endl;
			}
			if ((bestd->breakindex == int(bestd->seq.size()) - 1) && (ins->t[bestd->index].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "break pulled" << endl;
				pull_break(sol, bestd->index);
			}
			if ((beste->breakindex == int(beste->seq.size()) - 1) && (ins->t[beste->index].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "break pulled" << endl;
				pull_break(sol, beste->index);
			}
			//sol.check();
			succes = false;
		}//end if improvement
	}//end while improvement
	return succes;
}//end move_nb

bool Moves::swap2_nb(Sol& sol,int mode)//swap 2 vertices from two distinct tours in order to save traveltime
{
	bool improvement = true;
	bool succes = false;
	int iter = 0;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 0.0001;
		Sol::Tour *bestd=NULL;
		Sol::Tour *beste=NULL;
		int besti;
		int bestj;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			Sol::Tour* tourd = &sol.tours[d];
			for (int i = 0; i < tourd->seq.size() - 2; ++i)
			{
				Ins::Vertex* w = tourd->seq[i];
				Ins::Vertex* x = tourd->seq[i + 1];//vertex from tour d that will be swapped to tour e
				Ins::Vertex* y = tourd->seq[i + 2];
				double ttwxy = (tourd->deptime[i + 2] - (y->serv + tourd->action[i + 2] * ins->breakdur + x->serv + tourd->action[i + 1] * ins->breakdur)) - tourd->deptime[i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					Sol::Tour* toure = &sol.tours[e];
					if (d != e)//no relocate on the same tour
					{
						for (int j = 0; j < toure->seq.size() - 2; ++j)
						{
							Ins::Vertex* a = toure->seq[j];
							Ins::Vertex* b = toure->seq[j + 1];//vertex from tour e that will be swapped to tour d
							Ins::Vertex* c = toure->seq[j + 2];
							if ((a->nbi[e][x->index]) && (x->nbi[e][b->index]) && (w->nbi[d][y->index]))
							{
								if ((toure->weight + (x->weight-b->weight) <= ins->t[e].W_max) && (toure->volume + (x->weight-b->volume) <= ins->t[e].V_max)&&(tourd->weight + (b->weight - x->weight) <= ins->t[d].W_max)&&(tourd->volume + (b->volume - x->volume) <= ins->t[d].V_max))
								{//check if a potential increase in volume and weight is allowed on each tour
									
									//local evaluation on path d
									//the break can come too early due to potential traveltime decrease 
									bool reqbreakd = false;
									if (tourd->breakindex >= i + 1)//update break if it is currently position at x or after
									{
										reqbreakd = true;
									}
									double departuretime = ins->t[d].EDT + tourd->deptime[i];
									// w to b
									double arrivaltime = ins->arrival_time(w->con[b->index], departuretime);
									if ((reqbreakd) && (max(b->LTW[d] - ins->breakdur,arrivaltime) >= ins->breakstart))
									{
										if(max(b->LTW[d] - ins->breakdur, arrivaltime)>ins->breakend)
										{
											continue;
										}
										reqbreakd = false;
										arrivaltime += ins->breakdur;
										
									}
									if (arrivaltime < b->LTW[d])
									{
										arrivaltime = b->LTW[d];
									}
									if (arrivaltime > b->UTW[d])
									{
										continue;
									}
									arrivaltime += b->serv;
									departuretime = arrivaltime;
									//traveltime from b to y (can be on the end depot)
									arrivaltime = ins->arrival_time(b->con[y->index], departuretime);
									if ((reqbreakd) && ((max(y->LTW[d] - ins->breakdur,arrivaltime) >= ins->breakstart) || (y->index == ins->maxvertices - 1)))
									{
										reqbreakd = false;
										arrivaltime += ins->breakdur;
									}
									if (arrivaltime < y->LTW[d])
									{
										arrivaltime = y->LTW[d];
									}
									if (arrivaltime > y->UTW[d])
									{
										continue;
									}
									double arrivaltimey =arrivaltime+ y->serv;
									double diffd = (ins->t[d].EDT+tourd->deptime[i + 2])-arrivaltimey;
									//local evaluation on path e
									//a break can come too early due to potential traveltime decrease 
									bool reqbreake = false;
									if (toure->breakindex >= j + 1)
									{
										reqbreake = true;
									}
									departuretime = ins->t[e].EDT + toure->deptime[j];
									//traveltime a to x
									arrivaltime = ins->arrival_time(a->con[x->index], departuretime);
									if ((reqbreake) && (max(x->LTW[e] - ins->breakdur,arrivaltime) >= ins->breakstart))
									{
										if (max(x->LTW[e] - ins->breakdur, arrivaltime) > ins->breakend)
										{
											continue;
										}
										reqbreake = false;
										arrivaltime += ins->breakdur;
									}
									if (arrivaltime < x->LTW[e])
									{
										arrivaltime = x->LTW[e];
									}
									if (arrivaltime > x->UTW[e])
									{
										continue;
									}
									arrivaltime += x->serv;
									departuretime = arrivaltime;
									//traveltime from x to c (can be the end depot)
									arrivaltime = ins->arrival_time(x->con[c->index], departuretime);
									if ((reqbreake) && ((max(c->LTW[e] - ins->breakdur,arrivaltime) >= ins->breakstart) || (c->index == ins->maxvertices - 1)))
									{
										reqbreake = false;
										arrivaltime += ins->breakdur;
									}
									if (arrivaltime< c->LTW[e])
									{
										arrivaltime = c->LTW[e];
									}
									if (arrivaltime > c->UTW[e])
									{
										continue;
									}
									double arrivaltimec =arrivaltime + c->serv;
									double diffe = (ins->t[e].EDT+toure->deptime[j + 2]) - arrivaltimec;
									double localdecreasetotal = diffd + diffe;
									//local improvement check: check if potential increase is allowed and whether there is an overall travel time gain
									if ((-diffd<= tourd->max_shift[i + 2]) && (-diffe <= toure->max_shift[j + 2]) && (localdecreasetotal > bestdecrease))
									{
										//check enddepot time on path d
										Tour temptourd = *tourd;
										temptourd.replace_vertex(b,i+1);
										double globaldecreasetotal = tourd->deptime.back() - temptourd.deptime.back();
										//check enddepot time on path e
										Tour temptoure = *toure;
										temptoure.replace_vertex(x,j+1);
										globaldecreasetotal += toure->deptime.back() - temptoure.deptime.back();
										//global improvement check
										if (globaldecreasetotal > bestdecrease)
										{
											bestdecrease = globaldecreasetotal;
											beste = toure;
											bestj = j+1;
											bestd = tourd;
											besti = i+1;
											improvement = true;
											if (mode == 0)
											{
												goto swap2;
											}
										}
									}//end local check
								}//end if still capacity free
							}//end if neighbour
						}//end for j
					}//end for d!=e
				}//end for path e
			}//end for master sol i
		}// end for masterpath d
		if (improvement)
		{
		swap2:
			++iter;
			Sol remember = sol;
			Ins::Vertex *x = bestd->seq[besti];
			Ins::Vertex *b = beste->seq[bestj];
			sol.replace_vertex(*bestd,b,besti);
			sol.replace_vertex(*beste,x,bestj);
			sol.available[bestd->seq[besti]->index] = false;
			double actualdecrease = 0.0;
			for (int t = 0; t < (int) sol.tours.size(); ++t)
			{
				actualdecrease += remember.tours[t].deptime.back()-sol.tours[t].deptime.back();
			}
			sol.check();
			if (abs(bestdecrease - actualdecrease) > 0.01)
			{
				cout << "error swap2" << endl;
			}
			succes = true;
		}//end if improvement
	}//end while improvement
	return succes;
}//end swap2
