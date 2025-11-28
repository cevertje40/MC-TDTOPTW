#include "moves.h"

using namespace std;

bool Moves::insert_nb(Sol& sol, int mode)//insert vertex into a tour in order to increase the score
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestRatio = 0.0;
		int position = -1;
		Ins::Vertex* candidate = NULL;
		Sol::Tour* besttour = NULL;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			int t = sol.tourindex[d];
			Sol::Tour* tour= &sol.tours[t];
			/*
			auto& seq = tour->seq;
			auto& action = tour->action;
			auto& deptime = tour->deptime;
			auto& maxShift = tour->max_shift;

			const double EDT = ins->t[t].EDT;
			const double breakDur = ins->breakdur;
			const double Vmax = ins->t[t].V_max;
			const double Wmax = ins->t[t].W_max;
			*/
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
					if ((sol.available[y->index]) && (y->nbi[t][z->index]))//y moet buur van z zijn
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
								if (ratiocheck > bestRatio)//enkel op minimale increase checken
								{//update candidates
									improvement = true;
									bestRatio = ratiocheck;
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
			//Sol remember = sol;
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

bool Moves::exchange_nb(Sol& sol, int mode)//replaces a vertex of a tour with non-included vertex in the same position in order to increase the score
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
									double ratiocheck = double(y->score - z->score) / max(1.0, (y->weight - z->weight));
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
			int endh = (int)sol.tours[t].seq.size();
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
}

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
			int endh = (int)sol.tours[t].seq.size();
			for (int h = 1; h < endh - 1; ++h)// for all included regular vertices in sol
			{
				//for all existing regular member vertices: remove 1 and revaluate solution
				Sol::Tour tourrem = sol.tours[sol.tourindex[d]];
				Ins::Vertex* r = tourrem.seq[h];//to be removed vertex
				tourrem.remove_vertex(h);
				int endj = (int)tourrem.seq.size();
				for (int j = 0; j < endj - 1; ++j)// for positions in tourrem
				{
					Ins::Vertex* x1 = tourrem.seq[j];//predecessor y
					Ins::Vertex* z1 = tourrem.seq[j + 1];//successor y
					int breakz1 = tourrem.action[j + 1];
					int nb_size = (int)x1->nb[t].size();
					for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
					{
						Ins::Vertex* y1 = x1->nb[t][i];//potential insertion at position j
						if ((sol.available[y1->index]) && (y1->nbi[t][z1->index]))//availability & nb check
						{
							if ((tourrem.weight + y1->weight <= ins->t[t].W_max) && (tourrem.volume + y1->volume <= ins->t[t].V_max))//cap constraint check
							{
								//gather departure time
								double currenttime = tourrem.deptime[j] + ins->t[t].EDT;//service bij x zit hier al in
								//travel time from x1 to y1
								double arrivaltime = ins->arrival_time(x1->con[y1->index], currenttime);
								if (arrivaltime < y1->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
								{
									arrivaltime = y1->LTW[t];
								}
								if (arrivaltime > y1->UTW[t])
								{
									continue;//infeasible
								}
								arrivaltime += y1->serv;
								double arrivaltimey = arrivaltime;
								//travel time from y1 to z
								arrivaltime = ins->arrival_time(y1->con[z1->index], arrivaltime);
								if (arrivaltime + breakz1 * (ins->breakdur) < z1->LTW[t])
								{
									arrivaltime = z1->LTW[t] - (breakz1 * ins->breakdur);
								}
								arrivaltime += z1->serv + breakz1 * ins->breakdur;
								double shift = (arrivaltime - ins->t[t].EDT) - tourrem.deptime[j + 1];//increase in travel time
								if (shift <= tourrem.max_shift[j + 1])//check if first vertex can be inserted
								{
									//execute the first insertion
									Sol::Tour tourremins = tourrem;
									tourremins.insert_vertex(y1, j);
									int endk = (int)tourremins.seq.size();
									//check if second vertex can be inserted in this new tour
									for (int k = 0; k < endk - 1; ++k)// for all positions in tourremins
									{
										Ins::Vertex* x2 = tourremins.seq[k];//predecessor y2
										Ins::Vertex* z2 = tourremins.seq[k + 1];//successor y2
										int breakz2 = tourremins.action[k + 1];
										int nb_size = (int)x2->nb[t].size();
										for (int i = 0; i < nb_size - 1; ++i)//for all neighbours of the included vertex (non-enddepot)
										{
											Ins::Vertex* y2 = x2->nb[t][i];//potential insertion at position k
											if ((y1->score + y2->score) > r->score)
											{
												if ((sol.available[y2->index]) && (y2->nbi[t][z2->index]) && (y2 != y1))//availability & nb check & two insertions need to be different vertices
												{
													if ((tourremins.weight + y2->weight <= ins->t[t].W_max) && (tourremins.volume + y2->volume <= ins->t[t].V_max))//cap constraint check
													{
														//gather departure time
														double currenttime = tourremins.deptime[k] + ins->t[t].EDT;//service bij x zit hier al in
														//travel time from x to y2
														double arrivaltime = ins->arrival_time(x2->con[y2->index], currenttime);
														if (arrivaltime < y2->LTW[t])//break inserten kan niet dus break kan ltw niet dichter brengen
														{
															arrivaltime = y2->LTW[t];
														}
														if (arrivaltime > y2->UTW[t])
														{
															continue;//infeasible
														}
														arrivaltime += y2->serv;
														//travel time from y2 to z
														arrivaltime = ins->arrival_time(y2->con[z2->index], arrivaltime);
														if (arrivaltime + breakz2 * (ins->breakdur) < z2->LTW[t])
														{
															arrivaltime = z2->LTW[t] - (breakz2 * ins->breakdur);
														}
														arrivaltime += z2->serv + breakz2 * ins->breakdur;
														double shift = (arrivaltime - ins->t[t].EDT) - tourremins.deptime[k + 1];//increase in travel time
														if (shift <= tourremins.max_shift[k + 1])//check of het punt geinsert kan worden
														{
															improvement = true;
															double ratio = (y1->score + y2->score) - r->score;
															if (ratio > bestratio)
															{
																bestratio = ratio;
																besttour = &sol.tours[sol.tourindex[d]];
																best_pos_rem = h;
																best_pos_ins1 = j;
																best_pos_ins2 = k;//index after first insertion
																bestcandidate1 = y1;
																bestcandidate2 = y2;
																goto one_two_replace;
															}//end if bestratio improved
														}//end if insertion 2
													}//end cap check 2
												}//end availability check
											}//end if score improvement
										}//end nb check 2
									}//end for all positions in tourremins 2
								}//end if shift succes first insertion
							}//end cap check 1
						}//end availability & nb check 1
					}//end nb 1
				}//end for all positions in tourrem 1
			}//end for all included vertices
		}//end for all tours
		if (improvement)
		{
			one_two_replace:
			sol.remove_vertex(*besttour, best_pos_rem);
			sol.insert_vertex(*besttour,bestcandidate1, best_pos_ins1);
			sol.insert_vertex(*besttour,bestcandidate2, best_pos_ins2);
			//sol.check();
			//cout << "debug here" << endl;
		}
	}//end while improvement
	return succes;
}//end one_two_replace_nb

bool Moves::or_opt(Sol& sol, int mode)
{
	bool succes = false;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		bool improvement = true;
		Sol::Tour& tour = sol.tours[d];
		int size = int(tour.seq.size());
		while (improvement)
		{
			improvement = false;
			double bestdecrease = 1e-9;
			int bestfrom_start;
			int bestfrom_end;
			int bestposition;
			int bestbreakindex;
			vector<Ins::Vertex*> bestsubsequence;
			for (int from_start = 1; from_start < size-2; ++from_start)
			{
				for (int from_end = from_start + 1; from_end < size - 1; ++from_end)
				{
					for (int to = 1; to < size - 1; ++to)
					{
						if (to >= from_start && to <= from_end) continue;  // skip invalid insertions
						
						if (to == from_start || to == from_end + 1) continue;// skip insert in same position 


						// Step 1: Extract the subsequence from the route
						vector<Ins::Vertex*> subsequence(tour.seq.begin() + from_start, tour.seq.begin() + from_end + 1);
						Sol::Tour tourtry = tour;
						// Step 2: Remove the subsequence from the original position in the route
						tourtry.seq.erase(tourtry.seq.begin() + from_start, tourtry.seq.begin() + from_end + 1);
						// Step 3: Insert the subsequence at the new position in the route
						int position = to;
						if (position > from_start)
						{
							// If inserting later in the route, adjust the position since the original segment has been removed
							position -= (from_end - from_start);
						}
						tourtry.seq.insert(tourtry.seq.begin() + position, subsequence.begin(), subsequence.end());
						double departuretime = ins->t[tour.index].EDT;
						double arrivaltime=DBL_MAX;
						int breakindex=-1;
						bool reqbreak = true;
						for (int i = 1; i < size; ++i)
						{
							Ins::Vertex* o = tourtry.seq[i - 1];
							Ins::Vertex* p = tourtry.seq[i];
							arrivaltime = ins->arrival_time(o->con[p->index], departuretime);
							if ((reqbreak) && ((max(p->LTW[d] - ins->breakdur, arrivaltime) >= ins->breakstart) || (p->index == ins->maxvertices - 1)))
							{
								if (arrivaltime > ins->breakend)
								{
									arrivaltime = DBL_MAX;
									break;
								}
								breakindex = i;
								reqbreak = false;
								arrivaltime += ins->breakdur;
							}
							if (arrivaltime < p->LTW[d])
							{
								arrivaltime = p->LTW[d];
							}
							if (arrivaltime > p->UTW[d])
							{
								arrivaltime = DBL_MAX;
								break;
							}
							arrivaltime += p->serv;
							departuretime = arrivaltime;
						}
						arrivaltime -= ins->t[tour.index].EDT;
						double decrease = tour.deptime.back()-arrivaltime;
						if (decrease > bestdecrease)
						{
							improvement = true;
							bestdecrease = decrease;
							bestfrom_start = from_start;
							bestfrom_end = from_end;
							bestposition = position;
							bestsubsequence = subsequence;
							bestbreakindex = breakindex;
						}//end if
					}//end to
				}//end from end
			}//end from start
			if (improvement)
			{
				Sol remember = sol;
				sol.tours[d].seq.erase(sol.tours[d].seq.begin() + bestfrom_start, sol.tours[d].seq.begin() + bestfrom_end + 1);
				sol.tours[d].seq.insert(sol.tours[d].seq.begin() + bestposition, bestsubsequence.begin(), bestsubsequence.end());
				sol.tours[d].update_break(bestbreakindex);
				/*
				double actualdecrease= remember.tours[d].deptime.back() - sol.tours[d].deptime.back();
				if (!sol.check())
				{
					cout << "error in or-opt" << endl;
				}
				if (fabs(bestdecrease - actualdecrease) > 0.01)
				{
					cout << "error or-opt" << endl;
				}
				*/
				succes = true;
			}
		}
	}
	return succes;
}//end or_opt


static void dbg_dump_break_windows(const Ins* ins, const Tour& t, int d) {
	const int n = (int)t.seq.size();
	const double EDT = ins->t[d].EDT;
	const double LAT = ins->t[d].LAT;
	const double B = ins->breakdur;
	const double BEND = ins->breakend;
	const int ENDDEP = ins->maxvertices - 1;

	// Build arr0 (arrival-before-wait) on the CURRENT schedule clock
	std::vector<double> arr0(n, 0.0);
	{
		double cur = EDT;
		for (int u = 0; u < n - 1; ++u) {
			auto* o = t.seq[u];
			auto* p = t.seq[u + 1];
			double a = t.ins->arrival_time(o->con[p->index], cur);
			arr0[u + 1] = a;
			double s = (a < p->LTW[d] ? p->LTW[d] : a);
			cur = s + p->serv; // no-break semantics just for diagnostics
		}
	}

	std::cerr << "  B=" << B << " breakstart=" << ins->breakstart
		<< " breakend=" << BEND << " LAT=" << LAT << "\n";

	std::cerr << "  u | arr0 | earliest | latest | margin | LTW..UTW | isDepot\n";
	for (int u = 1; u < n; ++u) {
		const bool depot = (t.seq[u]->index == ENDDEP);
		double earliest = depot ? arr0[u] : std::max(arr0[u], ins->breakstart);
		double latest = depot ? std::min(BEND, LAT - B)
			: std::min(BEND, t.seq[u]->UTW[d] - B);
		double margin = latest - earliest;
		std::cerr << "  " << u
			<< " | " << arr0[u]
			<< " | " << earliest
			<< " | " << latest
			<< " | " << margin
			<< " | [" << t.seq[u]->LTW[d] << "," << t.seq[u]->UTW[d] << "]"
			<< " | " << depot << "\n";
	}
}

template<RatioFn RATIO>
boost::heap::priority_queue<One_one_rep_nb>
Moves::one_one_replace_gen_nb_kernel(Sol& sol, TabuVector& tabulist, int globalbest)
{
	boost::heap::priority_queue<One_one_rep_nb> adm_nb;

#ifdef _OPENMP
	const int max_threads = omp_get_max_threads();
#else
	const int max_threads = 1;
#endif
	std::vector<boost::heap::priority_queue<One_one_rep_nb>> tls_queues(max_threads);

#pragma omp parallel
	{
#ifdef _OPENMP
		const int tid = omp_get_thread_num();
#else
		const int tid = 0;
#endif
		auto& local_q = tls_queues[tid];

#pragma omp for nowait
		for (int d = 0; d < ins->maxtours; ++d)
		{
			double local_bestkey = -DBL_MAX;

			// ===== INSERT PART =====
			Sol::Tour& tour = sol.tours[d];
			const auto& Td = ins->t[d];
			const double EDT = Td.EDT;
			const double Wmax = Td.W_max, Vmax = Td.V_max, Tmax = Td.T_max;
			const double B = ins->breakdur;

			const int n = (int)tour.seq.size();
			for (int j = 0; j < n - 1; ++j)
			{
				Ins::Vertex* x = tour.seq[j];
				Ins::Vertex* z = tour.seq[j + 1];
				const int breakz = tour.action[j + 1];

				const double t0 = tour.deptime[j] + EDT;
				const double base_up = ins->arrival_time(x->con[z->index], t0) - t0;

				// projected break index after insertion at boundary j+1
				const int rb = tour.breakindex;          // in current tour
				const int nrem = n;                      // pre-insertion length
				const int r_post = (j + 1 <= rb) ? (rb + 1) : rb;  // can be == nrem
				const bool to_depot = (r_post == nrem);
				const int  rref = to_depot ? (nrem - 1) : r_post;   // safe index for arrays

				// pre-wait between (j+1) and the projected break position
				double pre_wait = 0.0;
				if (to_depot) pre_wait = tour.pref_wait.back() - tour.pref_wait[j + 1];
				else if (r_post > j + 1) pre_wait = tour.pref_wait[r_post] - tour.pref_wait[j + 1];

				// margin to legally START a break at the projected node
				const double br_margin_at = tour.br_margin[rref];

				// total break budget available from boundary j+1 point of view
				const double break_budget = pre_wait + br_margin_at;

				// penalty if the break slides earlier (i.e., r_post != rb)
				const double penalty = (r_post != rb) ? std::max(0.0, B - tour.wait_at[rref]) : 0.0;
				const double post_budget = std::max(0.0, tour.max_shift[j + 1] - penalty);

				int nb_size = (int)x->nb[d].size();
				for (int i = 0; i < nb_size - 1; ++i) // skip depot
				{
					Ins::Vertex* y = x->nb[d][i];

					if (sol.score + y->score <= globalbest + 1e-9)
						if (tabulist.isTabu(y->index, d)) continue;

					if (!(sol.available[y->index] && x->nbi[d][y->index] && y->nbi[d][z->index])) continue;
					if (tour.weight + y->weight > Wmax + 1e-9) continue;
					if (tour.volume + y->volume > Vmax + 1e-9) continue;

					const double ins_lb = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
					const double dt_lb = std::max(0.0, ins_lb - base_up);

					// quick key ub + suffix wait gate (suffix measured from j+1)
					const double suffix_wait = tour.pref_wait.back() - tour.pref_wait[j + 1];
					double key_ub = RATIO(dt_lb, y->score, y->weight, y->volume, Tmax, Wmax, Vmax);
					if (key_ub <= local_bestkey + 1e-9) continue;
					// cheap permissive screen (unchanged)
					if (dt_lb > suffix_wait + tour.max_shift[j + 1] + 1e-9) continue;

					// break-aware local gates at boundary j+1
					if (dt_lb > break_budget+ post_budget + 1e-9) continue;

					// exact timing
					double at = ins->arrival_time(x->con[y->index], t0);
					if (at < y->LTW[d]) at = y->LTW[d];
					if (at > y->UTW[d] + 1e-9) continue;
					at += y->serv;

					at = ins->arrival_time(y->con[z->index], at);
					if (breakz)
					{
						const bool endp = (z->index == ins->maxvertices - 1);
						const double latest_start = endp
							? std::min(ins->breakend, Td.LAT - B)
							: std::min(ins->breakend, z->UTW[d] - B);

						if (!endp && at < ins->breakstart) at = ins->breakstart;
						if (at > latest_start + 1e-9) continue; // break would start too late
						at += B;
					}

					if (at < z->LTW[d]) at = z->LTW[d];
					if (at > z->UTW[d] + 1e-9) continue;
					at += z->serv;


					bool survives_downstream = true;
					for (int p = j + 1; p < n - 1 && survives_downstream; ++p) {
						double Wprefix = tour.pref_wait[p + 1] - tour.pref_wait[j + 1];
						if (Wprefix < 0.0) Wprefix = 0.0;

						double residual = dt_lb - Wprefix;
						if (residual < 0.0) residual = 0.0;

						if (residual > tour.max_shift[p + 1] + 1e-9)
							survives_downstream = false;
					}
					// boundary j
					if (survives_downstream) {
						double baseUB_j = ins->arrival_time(x->con[y->index], t0) - t0;
						double insLB_j = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
						double dtLB_j = std::max(0.0, insLB_j - baseUB_j);

						if (dtLB_j > tour.max_shift[j + 1] + 1e-9)
							survives_downstream = false;
					}
					if (!survives_downstream) continue;

					const double shift = (at - EDT) - tour.deptime[j + 1];
					// same break-aware budgets apply to the exact shift
					if (shift > break_budget + post_budget + 1e-9) continue;


					// Compute global end-cap (absolute time) for this tour
					const double end_cap_abs = std::min({ ins->t[d].LAT,
														  ins->breakend + ins->breakdur,
														  ins->t[d].EDT + ins->t[d].T_max });

					// Current absolute route end (before applying this move) at end depot:
					const double cur_end_abs = ins->t[d].EDT + tour.deptime.back();

					// Remaining absolute slack now:
					const double global_slack_now = end_cap_abs - cur_end_abs;

					// After we apply this insertion locally, the route shifts by `shift`.
					// To still have room for *one* break of length B anywhere, we require:
					if (shift > global_slack_now - ins->breakdur + 1e-9) {
						continue; // reject this move globally
					}

					const double key = RATIO(shift, y->score, y->weight, y->volume, Tmax, Wmax, Vmax);
					if (key > local_bestkey + 1e-9)
					{
						local_bestkey = key;
						local_q.push(One_one_rep_nb(d, -1, j, y, y->score, key, { {}, { y } }));


						// --- DEBUG: verify break can be placed after first insertion ---
						{
							Sol::Tour t = tour;                 // for insert-only path; use tourrem in replace path
							t.insert_vertex(y, j);              // or y1 for one-two, etc.
							t.update_break();                   // recompute placement
							bool ok = t.check();                // your checker

							if (!ok || t.breakindex < 0) {
								std::cerr << "[FAIL] insert-only feasibility failed after update_break\n";
								std::cerr << "  d=" << d << " j=" << j
									<< " B=" << ins->breakdur
									<< " rb(before)=" << tour.breakindex
									<< " breakindex(after)=" << t.breakindex << "\n";

								auto dump_vec = [&](const char* name, const std::vector<double>& v) {
									std::cerr << "  " << name << ":";
									for (size_t ii = 0; ii < v.size(); ++ii) std::cerr << " " << v[ii];
									std::cerr << "\n";
									};
								auto dump_vec_i = [&](const char* name, const std::vector<int>& v) {
									std::cerr << "  " << name << ":";
									for (size_t ii = 0; ii < v.size(); ++ii) std::cerr << " " << v[ii];
									std::cerr << "\n";
									};

								dump_vec_i("action(before)", tour.action);
								dump_vec("max_shift(before)", tour.max_shift);
								dump_vec("br_margin(before)", tour.br_margin);
								dump_vec("wait_at(before)", tour.wait_at);
								dump_vec("pref_wait(before)", tour.pref_wait);

								dump_vec_i("action(after)", t.action);
								dump_vec("deptime(after)", t.deptime);
								dump_vec("max_shift(after)", t.max_shift);
								dump_vec("br_margin(after)", t.br_margin);
								dump_vec("wait_at(after)", t.wait_at);
								dump_vec("pref_wait(after)", t.pref_wait);

								// Also print UTW/LTW at/after the projected break area
								int N = (int)t.seq.size();
								std::cerr << "  LTW/UTW(after):";
								for (int ii = 0; ii < N; ++ii) {
									std::cerr << " [" << ii << ":" << t.seq[ii]->LTW[d] << "," << t.seq[ii]->UTW[d] << "]";
								}
								std::cerr << "\n";
							}
						}



					}
				}
			}

			const int endh = (int)sol.tours[d].seq.size();
			for (int h = 1; h < endh - 1; ++h)
			{
				Tour tourrem = sol.tours[d];
				Ins::Vertex* r = tourrem.seq[h];
				tourrem.remove_vertex(h);  // recomputes break & budgets (deptime/max_shift/pref_wait/br_margin/...)

				const auto& Td2 = ins->t[d];
				const double EDT2 = Td2.EDT;
				const double Wmax = Td2.W_max, Vmax = Td2.V_max, Tmax = Td2.T_max;
				const double B = ins->breakdur;
				const int    nrem = (int)tourrem.seq.size();
				const int    Nminus1 = nrem - 1;
				const double EPS = 1e-9;

				for (int j = 0; j < nrem - 1; ++j)
				{
					Ins::Vertex* x = tourrem.seq[j];
					Ins::Vertex* z = tourrem.seq[j + 1];
					const int breakz = tourrem.action[j + 1];

					const double t0 = tourrem.deptime[j] + EDT2;
					const double base_up = ins->arrival_time(x->con[z->index], t0) - t0;

					// projected break index after insertion at boundary j+1
					const int rb = tourrem.breakindex;             // current break index after removal
					const int r_post = (j + 1 <= rb) ? (rb + 1) : rb;  // may equal nrem
					const bool to_depot = (r_post == nrem);
					const int  rref = to_depot ? (nrem - 1) : r_post;

					// pre-wait between (j+1) and the projected break position
					double pre_wait = 0.0;
					if (to_depot) pre_wait = tourrem.pref_wait.back() - tourrem.pref_wait[j + 1];
					else if (r_post > j + 1) pre_wait = tourrem.pref_wait[r_post] - tourrem.pref_wait[j + 1];

					// margin to START a break at projected node
					const double br_margin_at = tourrem.br_margin[rref];

					// total break budget seen from boundary j+1
					const double break_budget = pre_wait + br_margin_at;

					// penalty if break slides earlier
					const double penalty = (r_post != rb) ? std::max(0.0, B - tourrem.wait_at[rref]) : 0.0;
					const double post_budget = std::max(0.0, tourrem.max_shift[j + 1] - penalty);

					int nb_size = (int)x->nb[d].size();
					for (int i = 0; i < nb_size - 1; ++i) // skip depot
					{
						Ins::Vertex* y = x->nb[d][i];

						// score/tabu
						if (sol.score + (y->score - r->score) <= globalbest + EPS)
							if (tabulist.isTabu(r->index, d) || tabulist.isTabu(y->index, d)) continue;

						// adjacency + capacities
						if (!(sol.available[y->index] && x->nbi[d][y->index] && y->nbi[d][z->index])) continue;
						if (tourrem.weight + y->weight > Wmax + EPS) continue;
						if (tourrem.volume + y->volume > Vmax + EPS) continue;

						// deterministic lower bound for insert
						const double ins_lb = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
						const double dt_lb = std::max(0.0, ins_lb - base_up);

						// cheap permissive key bound + suffix-wait screen
						const double suffix_wait = tourrem.pref_wait.back() - tourrem.pref_wait[j + 1];
						const double key_ub = RATIO(dt_lb, y->score - r->score, y->weight - r->weight, y->volume - r->volume,
							Tmax, Wmax, Vmax);
						if (key_ub <= local_bestkey + EPS) continue;
						if (dt_lb > suffix_wait + tourrem.max_shift[j + 1] + EPS) continue;

						// break-aware local budgets
						if (dt_lb > break_budget + post_budget + EPS) continue;

						// ----- exact timing through y and potential break at z -----
						double at = ins->arrival_time(x->con[y->index], t0);
						if (at < y->LTW[d]) at = y->LTW[d];
						if (at > y->UTW[d] + EPS) continue;
						at += y->serv;

						at = ins->arrival_time(y->con[z->index], at);
						if (breakz)
						{
							const bool endp = (z->index == ins->maxvertices - 1);
							const double latest_start = endp
								? std::min(ins->breakend, Td2.LAT - B)
								: std::min(ins->breakend, z->UTW[d] - B);

							if (!endp && at < ins->breakstart) at = ins->breakstart;
							if (at > latest_start + EPS) continue;
							at += B;
						}

						if (at < z->LTW[d]) at = z->LTW[d];
						if (at > z->UTW[d] + EPS) continue;
						at += z->serv;

						// exact local shift (after LTW snaps and potential break)
						const double shift = (at - EDT2) - tourrem.deptime[j + 1];

						// strong local break-aware budget on the *exact* shift
						if (shift > break_budget + post_budget + EPS) continue;

						// ---------- STRICT DOWNSTREAM SWEEP USING EXACT SHIFT ----------
						bool survives_downstream = true;

						// 1) boundary j itself: no prefix waiting before (j+1)
						if (shift > tourrem.max_shift[j + 1] + EPS)
							survives_downstream = false;

						// 2) boundaries after the insertion point
						for (int p = j + 1; p < Nminus1 && survives_downstream; ++p)
						{
							// waiting we can consume between (j+1) and (p+1)
							double Wprefix = tourrem.pref_wait[p + 1] - tourrem.pref_wait[j + 1];
							if (Wprefix < 0.0) Wprefix = 0.0;

							// residual delay reaching boundary (p+1)
							double residual = shift - Wprefix;
							if (residual < 0.0) residual = 0.0;

							// must fit entirely in local max_shift
							if (residual > tourrem.max_shift[p + 1] + EPS) {
								survives_downstream = false;
								break;
							}

							// UTW front-cap at vertex (p+1): new start-of-service must not exceed its local cap
							auto* vp1 = tourrem.seq[p + 1];
							double latest_svc_at_p1 = std::min(
								(double)vp1->UTW[d],
								EDT2 + tourrem.deptime[p + 1] + tourrem.max_shift[p + 1]
							);
							double new_svc_at_p1 = EDT2 + tourrem.deptime[p + 1] + residual;
							if (new_svc_at_p1 > latest_svc_at_p1 + EPS) {
								survives_downstream = false;
								break;
							}
						}
						if (!survives_downstream) continue;

						// --------- Global end-cap (absolute) ----------
						const double end_cap_abs = std::min({
							ins->t[d].LAT,
							ins->breakend + ins->breakdur,
							ins->t[d].EDT + ins->t[d].T_max
							});
						const double cur_end_abs = ins->t[d].EDT + tourrem.deptime.back();
						const double global_slack_now = end_cap_abs - cur_end_abs;
						if (shift > global_slack_now - B + EPS) continue;

						// key and enqueue
						const double key = RATIO(shift, y->score - r->score, y->weight - r->weight, y->volume - r->volume,
							Tmax, Wmax, Vmax);
						if (key <= local_bestkey + EPS) continue;

						local_bestkey = key;
						local_q.push(One_one_rep_nb(d, h, j, y, y->score - r->score, key, { { r }, { y } }));

						// --------- Optional DEBUG POST-CHECK (kept from your version) ---------
						{
							Sol::Tour dbg = sol.tours[d];
							dbg.remove_vertex(h);
							Sol::Tour old = dbg; // state after removal (for logging)

							// insert y after boundary j in the post-removal tour
							dbg.seq.insert(dbg.seq.begin() + j + 1, y);
							dbg.action.insert(dbg.action.begin() + j + 1, 0);
							dbg.deptime.insert(dbg.deptime.begin() + j + 1, 0);
							dbg.max_shift.insert(dbg.max_shift.begin() + j + 1, 0);

							dbg.score += y->score;
							dbg.weight += y->weight;
							dbg.volume += y->volume;

							bool ok = dbg.update_break();
							dbg.calc_maxshift();

							if (!dbg.check())
							{
								std::cerr << "[FAIL] replace feasibility failed after update_break\n";
								std::cerr << "  d=" << d << " j=" << j
									<< " B=" << ins->breakdur
									<< " rb(before)=" << old.breakindex
									<< " breakindex(after)=" << dbg.breakindex << "\n";

								auto dump_vec = [&](const char* name, const std::vector<double>& v) {
									std::cerr << "  " << name << ":";
									for (size_t ii = 0; ii < v.size(); ++ii) std::cerr << " " << v[ii];
									std::cerr << "\n";
									};
								auto dump_vec_i = [&](const char* name, const std::vector<int>& v) {
									std::cerr << "  " << name << ":";
									for (size_t ii = 0; ii < v.size(); ++ii) std::cerr << " " << v[ii];
									std::cerr << "\n";
									};

								dump_vec_i("action(before)", old.action);
								dump_vec("max_shift(before)", old.max_shift);
								dump_vec("br_margin(before)", old.br_margin);
								dump_vec("wait_at(before)", old.wait_at);
								dump_vec("pref_wait(before)", old.pref_wait);

								dump_vec_i("action(after)", dbg.action);
								dump_vec("deptime(after)", dbg.deptime);
								dump_vec("max_shift(after)", dbg.max_shift);
								dump_vec("br_margin(after)", dbg.br_margin);
								dump_vec("wait_at(after)", dbg.wait_at);
								dump_vec("pref_wait(after)", dbg.pref_wait);

								int N = (int)dbg.seq.size();
								std::cerr << "  LTW/UTW(after):";
								for (int ii = 0; ii < N; ++ii) {
									std::cerr << " [" << ii << ":" << dbg.seq[ii]->LTW[d]
										<< "," << dbg.seq[ii]->UTW[d] << "]";
								}
								std::cerr << "\n";
								dbg_dump_break_windows(ins, dbg, d);
							}
						}
					} // y loop
				} // j loop
			} // h loop
		} // for d
	} // parallel

	// Serial merge
	for (auto& q : tls_queues) {
		while (!q.empty()) {
			adm_nb.push(q.top());
			q.pop();
		}
	}
	return adm_nb;
}

boost::heap::priority_queue<One_one_rep_nb>Moves::one_one_replace_gen_nb(Sol& sol, TabuVector& tabulist, int globalbest, RatioKind kind)
{
	switch (kind) 
	{
		case RatioKind::SCORE:
			return one_one_replace_gen_nb_kernel<&ratio_scorediff>(sol, tabulist, globalbest);
		case RatioKind::TIME:
			return one_one_replace_gen_nb_kernel<&ratio_scorediff_time>(sol, tabulist, globalbest);
		case RatioKind::VOLUME:
			return one_one_replace_gen_nb_kernel<&ratio_scorediff_volume>(sol, tabulist, globalbest);
		case RatioKind::WEIGHT:
		default:
			return one_one_replace_gen_nb_kernel<&ratio_scorediff_weight>(sol, tabulist, globalbest);
	}
}

inline double base_travel_UB_over_interval(Ins* ins, Ins::Connec* c,double t0, double t1) 
{
	if (t1 < t0) { double tmp = t0; t0 = t1; t1 = tmp; }
	// Evaluate travel time at both ends and at every bucket boundary inside [t0, t1].
	double ub = 0.0;

	// convenience lambda-equivalent
	auto tt_at = [&](double t) -> double {
		return ins->arrival_time(c, t) - t; // UB at real departure t
		};

	ub = tt_at(t0);
	double v = tt_at(t1);
	if (v > ub) ub = v;

	int L = ins->find_t(t0);
	int R = ins->find_t(t1);
	// Scan interior boundaries: time_periods[q], q = L+1..R
	for (int q = L + 1; q <= R; ++q) {
		double tb = time_periods[q];                // boundary
		double val = tt_at(tb);                          // left-closed intervals → fine
		if (val > ub) ub = val;
	}
	return ub;
}





template<RatioFn RATIO>
boost::heap::priority_queue<Two_one_rep_nb>
Moves::two_one_replace_gen_nb_kernel(Sol& sol, TabuVector& tabulist, int globalbest)
{
	boost::heap::priority_queue<Two_one_rep_nb> adm_nb;

#ifdef _OPENMP
	const int max_threads = omp_get_max_threads();
#else
	const int max_threads = 1;
#endif
	std::vector<boost::heap::priority_queue<Two_one_rep_nb>> tls_queues(max_threads);

#pragma omp parallel
	{
#ifdef _OPENMP
		const int tid = omp_get_thread_num();
#else
		const int tid = 0;
#endif
		auto& local_q = tls_queues[tid];

#pragma omp for nowait
		for (int d = 0; d < ins->maxtours; ++d)
		{
			double local_bestkey = -DBL_MAX;

			const auto& Td = ins->t[d];
			const double EDT = Td.EDT;
			const double Wmax = Td.W_max, Vmax = Td.V_max, Tmax = Td.T_max;
			const double B = ins->breakdur;
			const double EPS = 1e-9;

			const int end = (int)sol.tours[d].seq.size();

			for (int g = 1; g < end - 1; ++g)            // first removed regular vertex
			{
				for (int h = g + 1; h < end - 1; ++h)    // second removed regular vertex
				{
					// Build tour with 2 removals (must recompute deptime/max_shift/pref_wait/br_margin/...).
					Sol::Tour tourrem = sol.tours[d];
					Ins::Vertex* r = tourrem.seq[g];
					Ins::Vertex* s = tourrem.seq[h];
					tourrem.remove_vertices(g, h);

					const int nrem = (int)tourrem.seq.size();
					const int Nminus1 = nrem - 1;

					for (int j = 0; j < nrem - 1; ++j) // insertion boundary in tourrem
					{
						Ins::Vertex* x = tourrem.seq[j];
						Ins::Vertex* z = tourrem.seq[j + 1];
						const int breakz = tourrem.action[j + 1];

						// Base timing at boundary
						const double t0 = tourrem.deptime[j] + EDT;
						const double base_up = ins->arrival_time(x->con[z->index], t0) - t0;

						// Project break index after insertion at boundary j+1
						const int rb = tourrem.breakindex; // current break in tourrem
						const int r_post_raw = (j + 1 <= rb) ? (rb + 1) : rb; // can be == nrem (to depot)
						const bool to_depot = (r_post_raw == nrem);
						const int  rref = to_depot ? (nrem - 1) : r_post_raw;

						// Waiting from boundary j+1 to the projected break
						double pre_wait = 0.0;
						if (to_depot) pre_wait = tourrem.pref_wait.back() - tourrem.pref_wait[j + 1];
						else if (r_post_raw > j + 1) pre_wait = tourrem.pref_wait[r_post_raw] - tourrem.pref_wait[j + 1];

						// Local break start margin at projected node
						const double br_margin_at = tourrem.br_margin[rref];

						// Combined “break budget” seen from j+1
						const double break_budget = pre_wait + br_margin_at;

						// Penalty if break slides earlier
						const double penalty = (r_post_raw != rb) ? std::max(0.0, B - tourrem.wait_at[rref]) : 0.0;
						const double post_budget = std::max(0.0, tourrem.max_shift[j + 1] - penalty);

						// Suffix wait from boundary j+1 (cheap gate)
						const double suffix_wait = tourrem.pref_wait.back() - tourrem.pref_wait[j + 1];

						// Neighbors
						const int nb_size = (int)x->nb[d].size();
						for (int i = 0; i < nb_size - 1; ++i) // skip end depot
						{
							Ins::Vertex* y = x->nb[d][i];

							// Tabu with aspiration vs globalbest
							if (sol.score + (y->score - (r->score + s->score)) <= globalbest + EPS)
							{
								if (tabulist.isTabu(r->index, d) || tabulist.isTabu(s->index, d) || tabulist.isTabu(y->index, d))
									continue;
							}

							// Availability & adjacency
							if (!(sol.available[y->index] && x->nbi[d][y->index] && y->nbi[d][z->index])) continue;

							// Capacities (after removals)
							if (tourrem.weight + y->weight > Wmax + EPS) continue;
							if (tourrem.volume + y->volume > Vmax + EPS) continue;

							// Deterministic LB gate (cheap)
							const double ins_lb = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
							const double dt_lb = std::max(0.0, ins_lb - base_up);

							// Cheap key ub & generic (suffix_wait + max_shift) gate
							const double key_ub = RATIO(
								dt_lb,
								y->score - (r->score + s->score),
								y->weight - (r->weight + s->weight),
								y->volume - (r->volume + s->volume),
								Tmax, Wmax, Vmax
							);
							if (key_ub <= local_bestkey + EPS) continue;
							if (dt_lb > suffix_wait + tourrem.max_shift[j + 1] + EPS) continue;

							// Break-aware local budgets at boundary j+1
							if (dt_lb > break_budget + post_budget + EPS) continue;

							// ----- Exact timing through y and potential break at z -----
							double at = ins->arrival_time(x->con[y->index], t0);
							if (at < y->LTW[d]) at = y->LTW[d];
							if (at > y->UTW[d] + EPS) continue;
							at += y->serv;

							at = ins->arrival_time(y->con[z->index], at);
							if (breakz)
							{
								const bool endp = (z->index == ins->maxvertices - 1);
								const double latest_start = endp
									? std::min(ins->breakend, Td.LAT - B)
									: std::min(ins->breakend, z->UTW[d] - B);

								if (!endp && at < ins->breakstart) at = ins->breakstart;
								if (at > latest_start + EPS) continue; // break too late
								at += B;
							}

							// LTW/UTW at z
							if (at < z->LTW[d]) at = z->LTW[d];
							if (at > z->UTW[d] + EPS) continue;
							at += z->serv;

							// Exact local shift (after snaps + possible break)
							const double shift = (at - EDT) - tourrem.deptime[j + 1];

							// Strong local break-aware budget on exact shift
							if (shift > break_budget + post_budget + EPS) continue;

							// ---------- STRICT DOWNSTREAM SWEEP USING EXACT SHIFT ----------
							bool survives_downstream = true;

							// 1) boundary j itself: no prefix waiting before (j+1)
							if (shift > tourrem.max_shift[j + 1] + EPS)
								survives_downstream = false;

							// 2) boundaries after insertion point
							for (int p = j + 1; p < Nminus1 && survives_downstream; ++p)
							{
								// waiting we can consume between (j+1) and (p+1)
								double Wprefix = tourrem.pref_wait[p + 1] - tourrem.pref_wait[j + 1];
								if (Wprefix < 0.0) Wprefix = 0.0;

								// residual delay reaching boundary (p+1)
								double residual = shift - Wprefix;
								if (residual < 0.0) residual = 0.0;

								// must fit entirely in local max_shift
								if (residual > tourrem.max_shift[p + 1] + EPS) {
									survives_downstream = false;
									break;
								}

								// UTW front-cap at vertex (p+1)
								auto* vp1 = tourrem.seq[p + 1];
								double latest_svc_at_p1 = std::min(
									(double)vp1->UTW[d],
									EDT + tourrem.deptime[p + 1] + tourrem.max_shift[p + 1]
								);
								double new_svc_at_p1 = EDT + tourrem.deptime[p + 1] + residual;
								if (new_svc_at_p1 > latest_svc_at_p1 + EPS) {
									survives_downstream = false;
									break;
								}
							}
							if (!survives_downstream) continue;

							// --------- Global end-cap (absolute) ----------
							const double end_cap_abs = std::min({
								ins->t[d].LAT,
								ins->breakend + ins->breakdur,
								ins->t[d].EDT + ins->t[d].T_max
								});
							const double cur_end_abs = ins->t[d].EDT + tourrem.deptime.back();
							const double global_slack_now = end_cap_abs - cur_end_abs;
							if (shift > global_slack_now - B + EPS) continue;

							// Key and enqueue
							const double key = RATIO(
								shift,
								y->score - (r->score + s->score),
								y->weight - (r->weight + s->weight),
								y->volume - (r->volume + s->volume),
								Tmax, Wmax, Vmax
							);
							if (key <= local_bestkey + EPS) continue;

							local_bestkey = key;
							local_q.push(Two_one_rep_nb(
								d, g, h, j, y,
								y->score - (r->score + s->score),
								key,
								{ { r, s }, { y } }
							));

							// ---------- Optional debug post-check ----------
							
							{
								Sol::Tour dbg = sol.tours[d];
								// remove s and r (in original order)
								dbg.remove_vertices(g, h);
								Sol::Tour old = dbg;

								// insert y after boundary j
								dbg.seq.insert(dbg.seq.begin() + j + 1, y);
								dbg.action.insert(dbg.action.begin() + j + 1, 0);
								dbg.deptime.insert(dbg.deptime.begin() + j + 1, 0);
								dbg.max_shift.insert(dbg.max_shift.begin() + j + 1, 0);

								dbg.score  += y->score;
								dbg.weight += y->weight;
								dbg.volume += y->volume;

								bool ok = dbg.update_break();
								dbg.calc_maxshift();
								if (!dbg.check()) {
									std::cerr << "[FAIL] two-one replace feasibility failed after update_break\n";
								}
							}
							
						} // neighbors
					} // j
				} // h
			} // g
		} // d
	} // omp parallel

	// Serial merge
	for (auto& q : tls_queues) {
		while (!q.empty()) {
			adm_nb.push(q.top());
			q.pop();
		}
	}
	return adm_nb;
}

boost::heap::priority_queue<Two_one_rep_nb>Moves::two_one_replace_gen_nb(Sol& sol, TabuVector& tabulist, int globalbest, RatioKind kind)
{
	switch (kind) 
	{
		case RatioKind::SCORE:
			return two_one_replace_gen_nb_kernel<&ratio_scorediff>(sol, tabulist, globalbest);
		case RatioKind::TIME:
			return two_one_replace_gen_nb_kernel<&ratio_scorediff_time>(sol, tabulist, globalbest);
		case RatioKind::VOLUME:
			return two_one_replace_gen_nb_kernel<&ratio_scorediff_volume>(sol, tabulist, globalbest);
		case RatioKind::WEIGHT:
		default:
			return two_one_replace_gen_nb_kernel<&ratio_scorediff_weight>(sol, tabulist, globalbest);
	}
}

template<RatioFn RATIO>
boost::heap::priority_queue<One_two_rep_nb>
Moves::one_two_replace_gen_nb_kernel(Sol& sol, TabuVector& tabulist, int globalbest)
{
	struct SchedPack {
		std::vector<double> deptime, max_shift, pref_wait, br_margin, wait_at;
		std::vector<int>    action;
		int breakindex = -1;
	};

	// Minimal forward simulation: insert y after boundary j into `base`,
	// place one legal break, and fill post-insert arrays into `out`.
	// (Same as your previous simulate, but no seq copy-out to keep it light.)
	static bool simulate_after_insert(
		const Sol::Tour & base, int j, Ins::Vertex * y, int d,
		const Ins * ins, double EDT, double LAT, double B, SchedPack & out)
	{
		const int n0 = (int)base.seq.size();
		const int n = n0 + 1;

		std::vector<Ins::Vertex*> seq(n);
		std::vector<int> action(n, 0);

		for (int i = 0; i <= j; ++i) { seq[i] = base.seq[i]; action[i] = 0; }
		seq[j + 1] = y; action[j + 1] = 0;
		for (int i = j + 1; i < n0; ++i) { seq[i + 1] = base.seq[i]; action[i + 1] = base.action[i]; }

		std::vector<double> deptime(n, 0.0), arr0(n, 0.0), svc(n, 0.0);
		double t = EDT;
		for (int u = 0; u < n - 1; ++u) {
			auto* o = seq[u];
			auto* p = seq[u + 1];
			double a = ins->arrival_time(o->con[p->index], t);
			arr0[u + 1] = a;
			double s = (a < p->LTW[d] ? p->LTW[d] : a);
			if (s > p->UTW[d] + 1e-9) return false;
			svc[u + 1] = s;
			double leave = s + p->serv;
			deptime[u + 1] = leave - EDT;
			t = leave;
		}

		std::vector<double> max_shift(n, 0.0), wait_at(n, 0.0),
			pref_wait(n, 0.0), br_margin(n, 0.0);

		// wait_at & pref_wait
		for (int u = 1; u < n; ++u) {
			const bool depot = (seq[u]->index == ins->maxvertices - 1);
			double anchor = depot ? arr0[u] : std::max(arr0[u], ins->breakstart);
			double usable = svc[u] - anchor; if (usable < 0.0) usable = 0.0;
			wait_at[u] = usable;
		}
		for (int u = 1; u < n; ++u) pref_wait[u] = pref_wait[u - 1] + wait_at[u];

		// br_margin
		for (int u = 1; u < n; ++u) {
			const bool depot = (seq[u]->index == ins->maxvertices - 1);
			double latest_start = depot ? std::min(ins->breakend, LAT - B)
				: std::min(ins->breakend, seq[u]->UTW[d] - B);
			double earliest_start = depot ? arr0[u] : std::max(arr0[u], ins->breakstart);
			double margin = latest_start - earliest_start; if (margin < 0.0) margin = 0.0;
			br_margin[u] = margin;
		}

		auto feasible_if_break_at = [&](int u)->bool {
			double cur = EDT;
			for (int s = 0; s < n - 1; ++s) {
				auto* o = seq[s]; auto* p = seq[s + 1];
				double arr = ins->arrival_time(o->con[p->index], cur);
				if (s + 1 == u) {
					const bool depot = (p->index == ins->maxvertices - 1);
					double latest_start = depot ? std::min(ins->breakend, LAT - B)
						: std::min(ins->breakend, p->UTW[d] - B);
					if (!depot && arr < ins->breakstart) arr = ins->breakstart;
					if (arr > latest_start + 1e-9) return false;
					arr += B;
				}
				if (arr < p->LTW[d]) arr = p->LTW[d];
				if (arr > p->UTW[d] + 1e-9) return false;
				arr += p->serv;
				cur = arr;
			}
			return true;
			};

		int breakindex = -1;
		for (int u = 1; u < n; ++u) {
			const bool depot = (seq[u]->index == ins->maxvertices - 1);
			double latest_start = depot ? std::min(ins->breakend, LAT - B)
				: std::min(ins->breakend, seq[u]->UTW[d] - B);
			double earliest_start = depot ? arr0[u] : std::max(arr0[u], ins->breakstart);
			if (earliest_start <= latest_start + 1e-9 && feasible_if_break_at(u)) {
				breakindex = u; action[u] = 1; break;
			}
		}
		if (breakindex < 0) return false;

		// Forward recompute with chosen break; also set a conservative max_shift
		{
			double cur = EDT;
			for (int s = 0; s < n - 1; ++s) {
				auto* o = seq[s]; auto* p = seq[s + 1];
				double arr = ins->arrival_time(o->con[p->index], cur);
				if (s + 1 == breakindex) {
					const bool depot = (p->index == ins->maxvertices - 1);
					if (!depot && arr < ins->breakstart) arr = ins->breakstart;
					arr += B;
				}
				if (arr < p->LTW[d]) arr = p->LTW[d];
				double earliest = arr;
				arr += p->serv;

				// recompute deptime
				deptime[s + 1] = arr - EDT;
				// conservative slack = UTW-bound on start-of-service
				double latest_svc = (double)p->UTW[d];
				double cur_svc = earliest;
				double slack = latest_svc - cur_svc; if (slack < 0.0) slack = 0.0;
				max_shift[s + 1] = slack;

				cur = arr;
			}
		}
		// End depot
		{
			double arr = EDT + deptime[n - 2];
			arr = ins->arrival_time(seq[n - 2]->con[seq[n - 1]->index], arr);
			if (action[n - 1]) {
				double latest_start_depot = std::min(ins->breakend, LAT - B);
				if (arr > latest_start_depot + 1e-9) return false;
				arr += B;
			}
			if (arr < seq[n - 1]->LTW[d]) arr = seq[n - 1]->LTW[d];
			if (arr > seq[n - 1]->UTW[d] + 1e-9) return false;
			arr += seq[n - 1]->serv;
			deptime[n - 1] = arr - EDT;
			// slack at depot (unused downstream, but set for completeness)
			max_shift[n - 1] = std::max(0.0, (double)seq[n - 1]->UTW[d] - (arr - seq[n - 1]->serv));
		}

		out.deptime = std::move(deptime);
		out.max_shift = std::move(max_shift);
		out.pref_wait = std::move(pref_wait);
		out.br_margin = std::move(br_margin);
		out.wait_at = std::move(wait_at);
		out.action = std::move(action);
		out.breakindex = breakindex;
		return true;
	}

	boost::heap::priority_queue<One_two_rep_nb> adm_nb;

#ifdef _OPENMP
	const int max_threads = omp_get_max_threads();
#else
	const int max_threads = 1;
#endif
	std::vector<boost::heap::priority_queue<One_two_rep_nb>> tls_queues(max_threads);

#pragma omp parallel
	{
#ifdef _OPENMP
		const int tid = omp_get_thread_num();
#else
		const int tid = 0;
#endif
		auto& local_q = tls_queues[tid];

#pragma omp for schedule(guided)
		for (int d = 0; d < ins->maxtours; ++d)
		{
			const auto& Td = ins->t[d];
			const double EDT = Td.EDT;
			const double Wmax = Td.W_max, Vmax = Td.V_max, Tmax = Td.T_max;
			const double B = ins->breakdur;
			const double EPS = 1e-9;

			double local_bestkey = -DBL_MAX;

			const int endh = (int)sol.tours[d].seq.size();
			for (int h = 1; h < endh - 1; ++h)
			{
				// remove r
				Sol::Tour tourrem = sol.tours[d];
				Ins::Vertex* r = tourrem.seq[h];
				tourrem.remove_vertex(h);

				const int nrem = (int)tourrem.seq.size();

				for (int j = 0; j < nrem - 1; ++j)
				{
					Ins::Vertex* x1 = tourrem.seq[j];
					Ins::Vertex* z1 = tourrem.seq[j + 1];
					const int breakz1 = tourrem.action[j + 1];

					const double t0_1 = EDT + tourrem.deptime[j];
					const double base_up1 = ins->arrival_time(x1->con[z1->index], t0_1) - t0_1;

					// project break wrt current tourrem
					const int rb = tourrem.breakindex;
					const int rpost = (j + 1 <= rb) ? (rb + 1) : rb;
					const bool to_depot = (rpost == nrem);
					const int  rref = to_depot ? (nrem - 1) : rpost;

					double pre_wait = 0.0;
					if (to_depot) pre_wait = tourrem.pref_wait.back() - tourrem.pref_wait[j + 1];
					else if (rpost > j + 1) pre_wait = tourrem.pref_wait[rpost] - tourrem.pref_wait[j + 1];

					const double br_margin_at = tourrem.br_margin[rref];
					const double before_budget = pre_wait + br_margin_at;
					const double penalty = (rpost != rb) ? std::max(0.0, B - tourrem.wait_at[rref]) : 0.0;
					const double after_budget = std::max(0.0, tourrem.max_shift[j + 1] - penalty);
					const double suffix_wait = tourrem.pref_wait.back() - tourrem.pref_wait[j + 1];

					const int nb1 = (int)x1->nb[d].size();
					for (int i1 = 0; i1 < nb1 - 1; ++i1)
					{
						Ins::Vertex* y1 = x1->nb[d][i1];

						if (sol.score + (y1->score - r->score) <= globalbest + EPS)
							if (tabulist.isTabu(r->index, d) || tabulist.isTabu(y1->index, d)) continue;

						if (!(sol.available[y1->index] && x1->nbi[d][y1->index] && y1->nbi[d][z1->index])) continue;
						if (tourrem.weight + y1->weight > Wmax + EPS) continue;
						if (tourrem.volume + y1->volume > Vmax + EPS) continue;

						const double insLB1 = x1->con[y1->index]->determin + y1->serv + y1->con[z1->index]->determin;
						const double dtLB1 = std::max(0.0, insLB1 - base_up1);

						double key_ub = RATIO(dtLB1, y1->score - r->score, y1->weight - r->weight, y1->volume - r->volume, Tmax, Wmax, Vmax);
						if (key_ub <= local_bestkey + EPS) continue;
						if (dtLB1 > suffix_wait + tourrem.max_shift[j + 1] + 1e-7) continue;
						if (dtLB1 > before_budget + after_budget + EPS) continue;

						// exact y1 timing
						double at1 = ins->arrival_time(x1->con[y1->index], t0_1);
						if (at1 < y1->LTW[d]) at1 = y1->LTW[d];
						if (at1 > y1->UTW[d] + EPS) continue;
						at1 += y1->serv;

						at1 = ins->arrival_time(y1->con[z1->index], at1);
						if (breakz1) {
							const bool endp = (z1->index == ins->maxvertices - 1);
							const double latest_start = endp ? std::min(ins->breakend, Td.LAT - B)
								: std::min(ins->breakend, z1->UTW[d] - B);
							if (!endp && at1 < ins->breakstart) at1 = ins->breakstart;
							if (at1 > latest_start + EPS) continue;
							at1 += B;
						}
						if (at1 < z1->LTW[d]) at1 = z1->LTW[d];
						if (at1 > z1->UTW[d] + EPS) continue;
						at1 += z1->serv;

						const double shift1 = (at1 - EDT) - tourrem.deptime[j + 1];
						if (shift1 > before_budget + after_budget + EPS) continue;

						// strict downstream sweep for shift1 on tourrem
						if (shift1 > tourrem.max_shift[j + 1] + EPS) continue;
						{
							bool ok = true;
							for (int p = j + 1; p < nrem - 1 && ok; ++p) {
								double Wprefix = tourrem.pref_wait[p + 1] - tourrem.pref_wait[j + 1];
								if (Wprefix < 0.0) Wprefix = 0.0;
								double residual1 = shift1 - Wprefix;
								if (residual1 < 0.0) residual1 = 0.0;

								if (residual1 > tourrem.max_shift[p + 1] + EPS) { ok = false; break; }

								auto* vp1 = tourrem.seq[p + 1];
								double latest_svc = std::min((double)vp1->UTW[d],
									EDT + tourrem.deptime[p + 1] + tourrem.max_shift[p + 1]);
								double new_svc = EDT + tourrem.deptime[p + 1] + residual1;
								if (new_svc > latest_svc + EPS) { ok = false; break; }
							}
							if (!ok) continue;
						}

						// post–y1 schedule
						SchedPack pack;
						if (!simulate_after_insert(tourrem, j, y1, d, ins, EDT, Td.LAT, B, pack)) {
							continue;
						}

						// helper to map indices to post–y1 sequence
						auto post_seq_at = [&](int i) -> Ins::Vertex* {
							if (i <= j)        return tourrem.seq[i];
							else if (i == j + 1) return y1;
							else               return tourrem.seq[i - 1];
							};

						const int n2 = (int)pack.deptime.size();
						for (int k = 0; k < n2 - 1; ++k)
						{
							Ins::Vertex* x2 = post_seq_at(k);
							Ins::Vertex* z2 = post_seq_at(k + 1);
							const int breakz2 = pack.action[k + 1];

							const double t0_2 = EDT + pack.deptime[k];
							const double base_up2 = ins->arrival_time(x2->con[z2->index], t0_2) - t0_2;

							// budgets at k+1 must use POST–y1 arrays (pack.*)
							const int rb2 = pack.breakindex;
							const int r2post = (k + 1 <= rb2) ? (rb2 + 1) : rb2;
							const bool to_dep2 = (r2post == n2);
							const int  rref2 = to_dep2 ? (n2 - 1) : r2post;

							double pre_wait2 = 0.0;
							if (to_dep2) pre_wait2 = pack.pref_wait.back() - pack.pref_wait[k + 1];
							else if (r2post > k + 1) pre_wait2 = pack.pref_wait[r2post] - pack.pref_wait[k + 1];

							const double br_margin2 = pack.br_margin[rref2];
							const double before_budget2 = pre_wait2 + br_margin2;

							const double penalty2 = (r2post != rb2) ? std::max(0.0, B - pack.wait_at[rref2]) : 0.0;
							const double after_budget2 = std::max(0.0, pack.max_shift[k + 1] - penalty2);

							const double suffix_wait2 = pack.pref_wait.back() - pack.pref_wait[k + 1];

							// neighbors for y2
							const int nb2 = (int)x2->nb[d].size();
							for (int i2 = 0; i2 < nb2 - 1; ++i2)
							{
								Ins::Vertex* y2 = x2->nb[d][i2];
								if (y2 == y1) continue;

								if (sol.score + ((y1->score + y2->score) - r->score) <= globalbest + EPS)
									if (tabulist.isTabu(r->index, d) || tabulist.isTabu(y1->index, d) || tabulist.isTabu(y2->index, d)) continue;

								if (!(sol.available[y2->index] && x2->nbi[d][y2->index] && y2->nbi[d][z2->index])) continue;

								// capacity after y1
								if (tourrem.weight + y1->weight + y2->weight > Wmax + EPS) continue;
								if (tourrem.volume + y1->volume + y2->volume > Vmax + EPS) continue;

								const double insLB2 = x2->con[y2->index]->determin + y2->serv + y2->con[z2->index]->determin;
								const double dtLB2 = std::max(0.0, insLB2 - base_up2);

								double key_ub2 = RATIO(shift1 + dtLB2,
									(y1->score + y2->score) - r->score,
									(y1->weight + y2->weight) - r->weight,
									(y1->volume + y2->volume) - r->volume,
									Tmax, Wmax, Vmax);
								if (key_ub2 <= local_bestkey + EPS) continue;
								if (dtLB2 > suffix_wait2 + pack.max_shift[k + 1] + 1e-7) continue;
								if (dtLB2 > before_budget2 + after_budget2 + EPS) continue;

								// exact timing for y2 on POST–y1 clock
								double at2 = ins->arrival_time(x2->con[y2->index], t0_2);
								if (at2 < y2->LTW[d]) at2 = y2->LTW[d];
								if (at2 > y2->UTW[d] + EPS) continue;
								at2 += y2->serv;

								at2 = ins->arrival_time(y2->con[z2->index], at2);
								if (breakz2) {
									const bool endp = (z2->index == ins->maxvertices - 1);
									const double latest_start2 = endp ? std::min(ins->breakend, Td.LAT - B)
										: std::min(ins->breakend, z2->UTW[d] - B);
									if (!endp && at2 < ins->breakstart) at2 = ins->breakstart;
									if (at2 > latest_start2 + EPS) continue;
									at2 += B;
								}
								if (at2 < z2->LTW[d]) at2 = z2->LTW[d];
								if (at2 > z2->UTW[d] + EPS) continue;
								at2 += z2->serv;

								// incremental delay of y2 vs POST–y1 schedule
								const double shift2 = (at2 - EDT) - pack.deptime[k + 1];
								if (shift2 > before_budget2 + after_budget2 + EPS) continue;

								// strict downstream sweep for shift2 using POST–y1 arrays
								if (shift2 > pack.max_shift[k + 1] + EPS) continue;
								{
									bool ok2 = true;
									for (int p = k + 1; p < n2 - 1 && ok2; ++p) {
										double Wpref = pack.pref_wait[p + 1] - pack.pref_wait[k + 1];
										if (Wpref < 0.0) Wpref = 0.0;
										double residual2 = shift2 - Wpref;
										if (residual2 < 0.0) residual2 = 0.0;

										if (residual2 > pack.max_shift[p + 1] + EPS) { ok2 = false; break; }

										auto* vp1 = post_seq_at(p + 1);
										double latest_svc = std::min((double)vp1->UTW[d],
											EDT + pack.deptime[p + 1] + pack.max_shift[p + 1]);
										double new_svc = EDT + pack.deptime[p + 1] + residual2;
										if (new_svc > latest_svc + EPS) { ok2 = false; break; }
									}
									if (!ok2) continue;
								}

								// Global end-cap based on POST–y1 end time
								{
									const double end_cap_abs = std::min({
										ins->t[d].LAT,
										ins->breakend + ins->breakdur,
										ins->t[d].EDT + ins->t[d].T_max
										});
									const double cur_end_abs = EDT + pack.deptime.back();
									const double slack_now = end_cap_abs - cur_end_abs;
									// conservative: added delay ≈ shift2 (shift1 already included in pack)
									if (shift2 > slack_now - B + EPS) continue;
								}

								const double key = RATIO(shift1 + shift2,
									(y1->score + y2->score) - r->score,
									(y1->weight + y2->weight) - r->weight,
									(y1->volume + y2->volume) - r->volume,
									Tmax, Wmax, Vmax);
								if (key <= local_bestkey + EPS) continue;

								local_bestkey = key;
								local_q.push(One_two_rep_nb(
									d, h, j, k, y1, y2,
									(y1->score + y2->score) - r->score,
									key,
									{ { r }, { y1, y2 } }
								));
							} // i2
						} // k
					} // i1
				} // j
			} // h
		} // d
	} // omp parallel

	for (auto& q : tls_queues) { while (!q.empty()) { adm_nb.push(q.top()); q.pop(); } }
	return adm_nb;
}

boost::heap::priority_queue<One_two_rep_nb>Moves::one_two_replace_gen_nb(Sol& sol, TabuVector& tabulist, int globalbest, RatioKind kind)
{
	switch (kind) 
	{
		case RatioKind::SCORE:
			return one_two_replace_gen_nb_kernel<&ratio_scorediff>(sol, tabulist, globalbest);
		case RatioKind::TIME:
			return one_two_replace_gen_nb_kernel<&ratio_scorediff_time>(sol, tabulist, globalbest);
		case RatioKind::VOLUME:
			return one_two_replace_gen_nb_kernel<&ratio_scorediff_volume>(sol, tabulist, globalbest);
		case RatioKind::WEIGHT:
		default:
			return one_two_replace_gen_nb_kernel<&ratio_scorediff_weight>(sol, tabulist, globalbest);
	}
}

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
	auto try_take_break = [&](Ins::Vertex* v, int droute, double& arr) -> bool 
		{
		const int END_DEPOT = ins->maxvertices - 1;
		const double B = ins->breakdur;
		const bool endp = (v->index == END_DEPOT);
		double start = endp ? arr : std::max(arr, ins->breakstart);
		if (start > ins->breakend + 1e-9) return false;
		double latest_start = std::min(ins->breakend, ins->t[droute].LAT - B);
		if (start > latest_start + 1e-9) return false;
		arr = start + B;
		return true;
		};

	bool succes = false;
	for (int d = 0; d < ins->maxtours; ++d)
	{
		Sol::Tour& tour = sol.tours[d];
		bool improvement = true;
		while (improvement)
		{
			double bestdelta = 0.0;
			int besti=-1;
			int bestj=-1;
			int bestbreakindex = tour.breakindex;
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

						double lb = k->con[z->index]->determin + z->con[l->index]->determin +m->con[y->index]->determin + y->con[n->index]->determin- (k->con[y->index]->determin + y->con[l->index]->determin +m->con[z->index]->determin + z->con[n->index]->determin);
						if (-lb >= + 1e-9)
						{
							double delta_tt = 0.0;
							double departuretime = ins->t[d].EDT + tour.deptime[i - 1];//bij punt voor y
							double currenttime = departuretime;
							double newtraveltime;
							double oldtraveltime;
							//calculate swapped indices
							auto map_idx = [&](int u)->int
								{
									if (u == i) return j;
									if (u == j) return i;
									return u;
								};
							for (int l = i - 1; l < j + 1; ++l)
							{
								Ins::Vertex* first = tour.seq[map_idx(l)];
								Ins::Vertex* second = tour.seq[map_idx(l + 1)];
								double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
								//do not account for breaks in local evaluation
								if (arrivaltime < second->LTW[d])
								{
									arrivaltime = second->LTW[d];
								}
								if (arrivaltime > second->UTW[d] + 1e-9)
								{
									currenttime = DBL_MAX;
									break;
								}
								arrivaltime += second->serv;
								currenttime = arrivaltime;
							}// end l
							newtraveltime = currenttime - departuretime;
							oldtraveltime = (ins->t[d].EDT + tour.deptime[j + 1]) - departuretime;
							delta_tt = oldtraveltime - newtraveltime;
							if (delta_tt > bestdelta+ 1e-9)//local evaluation
							{
								bool reqbreak = true;
								int breakindex = tour.breakindex;
								currenttime = ins->t[d].EDT;
								for (int l = 0; l < end - 1; ++l)//global evaluation
								{
									Ins::Vertex* first = tour.seq[map_idx(l)];
									Ins::Vertex* second = tour.seq[map_idx(l + 1)];//can be the end depot
									double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
									//account for break
									if (reqbreak)
									{
										if (try_take_break(second, d, arrivaltime)) 
										{
											reqbreak = false;
											breakindex = l + 1;
										}
									}
									if (arrivaltime < second->LTW[d])
									{
										arrivaltime = second->LTW[d];
									}
									if (arrivaltime > second->UTW[d] + 1e-9)
									{
										currenttime = DBL_MAX;
										break;
									}
									arrivaltime += second->serv;
									currenttime = arrivaltime;
								}// end l
								// If break still pending at the depot, combination is infeasible
								if (reqbreak) 
								{
									currenttime = DBL_MAX;
								}
								newtraveltime = currenttime - departuretime;
								oldtraveltime = ((ins->t[d].EDT + tour.deptime[end - 1]) - departuretime);
								delta_tt = oldtraveltime-newtraveltime;
								if (delta_tt >= bestdelta + 1e-9)//global evaluation
								{
									improvement = true;
									bestdelta = delta_tt;
									bestbreakindex = breakindex;
									besti = i;
									bestj = j;
									if (mode == 0)
									{
										goto swap;
									}
								}//end global evaluation
							}// end local evaluation
						}//end tid filter
					}//end neighbourhood check
				}//end j
			}// end i
			if (improvement)
			{
				swap:
				Sol remember = sol;
				tour.swap_vertices(besti,bestj,bestbreakindex);
				
				double actualdecrease = remember.tours[d].deptime.back() - tour.deptime.back();
				if (!sol.check())
				{
					cout<< term::fg(term::Color::red) << "sol error in swap" << endl;
				}
				if (fabs(bestdelta - actualdecrease) > 0.01)
				{
					cout<< term::fg(term::Color::red) << "time decrease error swap: "<< bestdelta - actualdecrease << endl;
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
			double bestdelta = 1e-9;
			int besti = -1;
			int bestj = -1;
			int bestbreakindex = -1;
			improvement = false;
			int end = (int)tour.seq.size();
			for (int i = 1; i < end - 1; ++i)// depots don't count
			{
				for (int j = i + 2; j < end - 1; ++j)// true 2-opt reversals only (segment length >= 3), the rest is taken care of by swap
				{
					Ins::Vertex* k = tour.seq[i - 1];//pred partner 1
					Ins::Vertex* y = tour.seq[i];//partner 1
					Ins::Vertex* l = tour.seq[i + 1];//suc partner 1
					Ins::Vertex* m = tour.seq[j - 1];//pred  partner 2
					Ins::Vertex* z = tour.seq[j];// partner 2
					Ins::Vertex* n = tour.seq[j + 1];//suc partner 2
					if ((k->nbi[d][z->index]) && (z->nbi[d][m->index]) && (l->nbi[d][y->index]) && (y->nbi[d][n->index]))
					{
						double lb = (k->con[y->index]->determin + z->con[n->index]->determin)- (k->con[z->index]->determin + y->con[n->index]->determin);
						if (lb > 1e-9)
						{
							//calculate reversed subpath
							auto map_idx = [&](int u)->int
								{
									if (u < i) return u;
									if (u > j) return u;
									return i + j - u;
								};

							bool reqbreak = false;
							if (i <= tour.breakindex)
							{//break is positioned after i, so might need to be repositioned
								reqbreak = true;
							}
							double delta_tt = 0;
							double departuretime = ins->t[d].EDT + tour.deptime[i - 1];
							double currenttime = departuretime;
							double arrivaltime;
							int breakindex = tour.breakindex;
							for (int l = i - 1; l < j + 1; ++l)
							{
								Ins::Vertex* first = tour.seq[map_idx(l)];
								Ins::Vertex* second = tour.seq[map_idx(l + 1)];
								arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
								if ((reqbreak) && ((max(second->LTW[d] - ins->breakdur, arrivaltime) >= ins->breakstart + 1e-9) || (second->index == ins->maxvertices - 1)))
								{
									if (reqbreak && max(second->LTW[d] - ins->breakdur, arrivaltime) > ins->breakend+ 1e-9)
									{
										currenttime = DBL_MAX;
										break;
									}
									arrivaltime += ins->breakdur;
									reqbreak = false;
									breakindex = l + 1;//record breakindex
								}
								if (arrivaltime < second->LTW[d])
								{
									arrivaltime = second->LTW[d];
								}
								if (arrivaltime > second->UTW[d] + 1e-9)
								{
									currenttime = DBL_MAX;
									break;
								}
								currenttime = arrivaltime + second->serv;
							}// end l
							if (reqbreak == true)//als je aankomt bij het einddepot en nog steeds geen break genomen hebt
							{
								currenttime += ins->breakdur;//breaktime bijtellen bij aankomst tijd bij einddepot
								breakindex = end - 1;//record breakindex
							}
							if (currenttime == DBL_MAX) continue;
							delta_tt = (ins->t[d].EDT + tour.deptime[j + 1]) - currenttime;
							if (delta_tt > +1e-9)//local evaluation
							{
								Tour temptour = tour;
								temptour.opt_vertices(i, j,breakindex);
								delta_tt = tour.deptime.back() - temptour.deptime.back();
								if (delta_tt > bestdelta + 1e-9)//global evaluation
								{
									improvement = true;
									bestdelta = delta_tt;
									bestbreakindex = breakindex;
									besti = i;
									bestj = j;
								}
							}//end executed 2 opt
						}//end neighborhood check
					}//end tid filter
				}// end for al j
			}// end for all i
			if (improvement)
			{
				Sol remember = sol;
				tour.opt_vertices(besti,bestj,bestbreakindex);
				double actualdecrease = remember.tours[d].deptime.back() - tour.deptime.back();
				
				if (!sol.check())
				{
					cout<< term::fg(term::Color::red) << "error sol two-opt" << endl;
				}
				if (fabs(bestdelta - actualdecrease) > 0.01)
				{
					cout<< term::fg(term::Color::red) << "error time decrase two-opt:"<<bestdelta-actualdecrease << endl;
				}
				
				succes = true;
			}
		}//end while improvement
	}//end for all paths
	return succes;
}

bool Moves::shift_nb(Sol& sol, int mode)
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 1e-9;
		Sol::Tour* bestd = NULL;
		int besti=-1;
		int bestj=-1;
		int bestbreakindex=-1;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			Sol::Tour* tourd = &sol.tours[d];
			int end = (int)tourd->seq.size();
			for (int i = 1; i < end - 2; ++i)
			{
				for (int j = i+1; j < end - 1; ++j)
				{
					// --- index mapping: view of sequence after moving i after j ---
					auto at = [&](int v) -> Ins::Vertex* 
						{
						if (v < i)             return tourd->seq[v];
						if (v >= i && v < j)   return tourd->seq[v + 1];
						if (v == j)             return tourd->seq[i];
						/* v > j */             return tourd->seq[v];
						};

					bool reqbreak = (tourd->breakindex >= i);
					int breakindex = tourd->breakindex;
					double departuretime = ins->t[d].EDT + tourd->deptime[i - 1];
					double arrivaltime=0.0;
					for (int v = i; v <= j + 1; ++v)
					{
						Ins::Vertex* o = at(v - 1);
						Ins::Vertex* p = at(v);
						arrivaltime = ins->arrival_time(o->con[p->index], departuretime);
						if ((reqbreak) && ((max(p->LTW[d] - ins->breakdur, arrivaltime) >= ins->breakstart) || (p->index == ins->maxvertices - 1)))
						{
							if (arrivaltime > ins->breakend)
							{
								arrivaltime = DBL_MAX;
								break;
							}
							reqbreak = false;
							breakindex = v;
							arrivaltime += ins->breakdur;
						}
						if (arrivaltime < p->LTW[d])
						{
							arrivaltime = p->LTW[d];
						}
						if (arrivaltime > p->UTW[d])
						{
							arrivaltime = DBL_MAX;
							break;
						}
						arrivaltime += p->serv;
						departuretime = arrivaltime;
					}
					double localdecreasetotal = tourd->deptime[j + 1] - (arrivaltime - ins->t[d].EDT);
					if (localdecreasetotal > bestdecrease)
					{
						departuretime = arrivaltime;
						for (int v = j+2; v < end ; ++v)
						{
							Ins::Vertex* o = at(v - 1);
							Ins::Vertex* p = at(v);
							arrivaltime = ins->arrival_time(o->con[p->index], departuretime);
							if ((reqbreak) && ((max(p->LTW[d] - ins->breakdur, arrivaltime) >= ins->breakstart) || (p->index == ins->maxvertices - 1)))
							{
								if (arrivaltime > ins->breakend)
								{
									arrivaltime = DBL_MAX;
									break;
								}
								breakindex = v;
								reqbreak = false;
								arrivaltime += ins->breakdur;
							}
							if (arrivaltime < p->LTW[d])
							{
								arrivaltime = p->LTW[d];
							}
							if (arrivaltime > p->UTW[d])
							{
								arrivaltime = DBL_MAX;
								break;
							}
							arrivaltime += p->serv;
							departuretime = arrivaltime;
						}
						double globaldecrease = tourd->deptime.back() - (arrivaltime - ins->t[d].EDT);
						if (globaldecrease > bestdecrease)
						{
							bestdecrease = globaldecrease;
							bestd = tourd;
							besti = i;
							bestj = j;
							bestbreakindex=breakindex;
							improvement = true;
							if (mode == 0)
							{
								goto shift;
							}
						}//end global check
					}//end local check
				}//end for j
			}//end for i
		}//end for all tours
		if (improvement)
		{
		shift:
			//Sol remember = sol;
			Ins::Vertex* moved = bestd->seq[besti];
			bestd->seq.erase(bestd->seq.begin() + besti);
			bestd->seq.insert(bestd->seq.begin() + bestj, moved);
			bestd->update_break(bestbreakindex);
			/*
			double actualdecrease = 0.0;
			for (int t = 0; t < (int)sol.tours.size(); ++t)
			{
				actualdecrease += remember.tours[t].deptime.back() - sol.tours[t].deptime.back();
			}
			if (fabs(bestdecrease - actualdecrease) > 0.01)
			{
				cout << "error shift" << endl;
			}
			*/
			if ((bestd->breakindex == int(bestd->seq.size()) - 1) && (ins->t[bestd->index].LAT > ins->breakend + ins->breakdur))
			{
				//cout << "break pulled" << endl;
				pull_break(sol, bestd->index);
			}
			//if (!sol.check())
			//{
				//cout << "error in shift_nb" << endl;
			//}
			succes = true;
		}//end if improvement
	}//end while improvement
	return succes;
}

bool Moves::relocate_nb(Sol& sol,int mode)//move vertex x from tour d to tour e in order to save travel time
{
	bool improvement = true;
	bool succes = false;
	auto try_take_break = [&](Ins::Vertex* v, int droute, double& arr) -> bool
	{
			const int END_DEPOT = ins->maxvertices - 1;
			const double B = ins->breakdur;
			const bool endp = (v->index == END_DEPOT);
			double start = endp ? arr : std::max(arr, ins->breakstart);
			if (start > ins->breakend + 1e-9) return false;
			double latest_start = std::min(ins->breakend, ins->t[droute].LAT - B);
			if (start > latest_start + 1e-9) return false;
			arr = start + B;
			return true;
	};
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 1e-9;
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
				double sos_y = tourd->deptime[i + 2] - (y->serv + tourd->action[i + 2] * ins->breakdur);
				double sos_y_no_xsrv = sos_y - (tourd->action[i + 1] * ins->breakdur + x->serv);
				double ttwxy = sos_y_no_xsrv - tourd->deptime[i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					if (d != e)
					{
						Sol::Tour* toure = &sol.tours[e];
						int r = toure->breakindex;
						for (int j = 0; j < int(toure->seq.size()) - 2; ++j)
						{
							int r_post = (j + 1 <= r) ? (r + 1) : r; // new break position after inserting x between a and b
							Ins::Vertex* a = toure->seq[j];
							Ins::Vertex* b = toure->seq[j + 1];
							if ((a->nbi[e][x->index]) && (x->nbi[e][b->index]) && (w->nbi[d][y->index]))
							{
								if ((toure->weight + x->weight <= ins->t[e].W_max + 1e-9) && (toure->volume + x->volume <= ins->t[e].V_max + 1e-9))
								{
									double save_d = w->con[x->index]->determin + x->con[y->index]->determin - w->con[y->index]->determin;         // removal gain on d
									double cost_e = a->con[x->index]->determin + x->con[b->index]->determin - a->con[b->index]->determin;         // insertion cost on e
									double lb = save_d - cost_e;                          // >0 is promising
									if (lb > 1e-9)
									{
										//local evaluation on path e, break remains unchanged
										double depA = ins->t[e].EDT + toure->deptime[j];
										double ttab = (toure->deptime[j + 1] - (b->serv + toure->action[j + 1] * ins->breakdur)) - toure->deptime[j];
										// a -> x
										double arr = ins->arrival_time(a->con[x->index], depA);
										if (arr < x->LTW[e]) arr = x->LTW[e];
										if (arr > x->UTW[e] + 1e-9) continue;
										arr += x->serv;

										// x -> b
										double arr_b = ins->arrival_time(x->con[b->index], arr);

										// break-at-b if r_post == j+2
										if (r_post == j + 2) 
										{
											if (!try_take_break(b, e, arr_b)) continue;
										}

										// LTW at b (after any break)
										if (arr_b < b->LTW[e]) arr_b = b->LTW[e];
										// UTW guard (cheap, defensive; max_shift bound *should* imply it)
										if (arr_b > b->UTW[e] + 1e-9) continue;

										// This is start-of-service at b in the candidate schedule
										double ttaxb_bstart = arr_b - depA;
										
										// And time after serving b (for carrying on)
										double arrivalb_svc_done = arr_b + b->serv;
										
										//local evaluation on path d
										//due to potential traveltime decrease break can come to early
										bool reqbreak = (tourd->breakindex >= i + 1);
										//calculate wy
										double depW = ins->t[d].EDT + tourd->deptime[i];

										// w -> y
										double arr_y = ins->arrival_time(w->con[y->index], depW);
										// if break moved earlier, try to take it at y (start-only semantics)
										if (reqbreak) 
										{
											if (try_take_break(y, d, arr_y)) 
											{
												reqbreak=false;
											}
										}

										// LTW at y
										if (arr_y < y->LTW[d]) arr_y = y->LTW[d];
										if (arr_y > y->UTW[d] + 1e-9) continue;

										// start-of-service at y (candidate)
										double ttwy_bstart = arr_y - depW;
										//double arrivaly = arrivaltime + y->serv;
										
										double shift = (arrivalb_svc_done - ins->t[e].EDT) - toure->deptime[j + 1];
										double localdecreasetotal = ttwxy + ttab - (ttwy_bstart + ttaxb_bstart);
										if ((shift <= toure->max_shift[j + 1] +  + 1e-9) && (localdecreasetotal > bestdecrease + 1e-9))//local improvement check
										{
											//check enddepot time on path d
											Tour temptourd = *tourd;
											temptourd.remove_vertex(i + 1);
											double globaldecreasetotal = (tourd->deptime.back() - temptourd.deptime.back());

											//check enddepot time on path e
											double currenttime = arrivalb_svc_done;
											for (int m = j + 2; m < (int)toure->seq.size(); ++m) 
											{
												Ins::Vertex* o = toure->seq[m - 1];
												Ins::Vertex* p = toure->seq[m];

												double arrm = ins->arrival_time(o->con[p->index], currenttime);

												// new index for p after insertion is (m + 1)
												int new_idx_p = m + 1;
												if (new_idx_p == r_post) 
												{
													if (!try_take_break(p, e, arrm)) 
													{
														currenttime = DBL_MAX; break; 
													}
												}

												if (arrm < p->LTW[e])
													arrm = p->LTW[e];
												

												arrm += p->serv;
												currenttime = arrm;
											}
											if (currenttime == DBL_MAX) continue;
											//cout << "enddepot time e: " << currenttime - t[e].EDT << endl;
											//cout << "increase e: " << (currenttime - t[e].EDT) - sol.traveltime[e].back() << endl;
											globaldecreasetotal -= ((currenttime - ins->t[e].EDT) - toure->deptime.back());
											if (globaldecreasetotal > bestdecrease + 1e-9)
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
								}//end tid check
							}//end if neighbour
						}//end for j
					}//end for d!=e
				}//end for tour e
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
			if (fabs(bestdecrease - actualdecrease) > 0.01)
			{
				cout<< term::fg(term::Color::red) << "error time decrease relocate"<<bestdecrease-actualdecrease << endl;
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
			if (!sol.check())
			{
				cout<< term::fg(term::Color::red) << "error sol in relocate" << endl;
			}
			succes = true;
		}//end if improvement
	}//end while improvement
	return succes;
}//end move_nb

bool Moves::swap2_nb(Sol& sol,int mode)//swap 2 vertices from two distinct tours in order to save traveltime
{
	bool improvement = true;
	bool succes = false;
	while (improvement)
	{
		improvement = false;
		double bestdecrease = 1e-9;
		Sol::Tour *bestd=NULL;
		Sol::Tour *beste=NULL;
		int bestbreakindexd = -1;
		int bestbreakindexe = -1;
		int besti=-1;
		int bestj=-1;
		for (int d = 0; d < ins->maxtours; ++d)
		{
			Sol::Tour* tourd = &sol.tours[d];
			for (int i = 0; i < tourd->seq.size() - 2; ++i)
			{
				Ins::Vertex* w = tourd->seq[i];
				Ins::Vertex* x = tourd->seq[i + 1];//vertex from tour d that will be swapped to tour e
				Ins::Vertex* y = tourd->seq[i + 2];
				for (int e = d+1; e < ins->maxtours; ++e)
				{
					Sol::Tour* toure = &sol.tours[e];
					for (int j = 0; j < toure->seq.size() - 2; ++j)
					{
						Ins::Vertex* a = toure->seq[j];
						Ins::Vertex* b = toure->seq[j + 1];//vertex from tour e that will be swapped to tour d
						Ins::Vertex* c = toure->seq[j + 2];
						if ((a->nbi[e][x->index]) && (x->nbi[e][c->index]) && (w->nbi[d][b->index])&&(b->nbi[d][y->index]))
						{
							if ((toure->weight + (x->weight-b->weight) <= ins->t[e].W_max + 1e-9) && (toure->volume + (x->volume-b->volume) <= ins->t[e].V_max + 1e-9)&&(tourd->weight + (b->weight - x->weight) <= ins->t[d].W_max + 1e-9)&&(tourd->volume + (b->volume - x->volume) <= ins->t[d].V_max + 1e-9))
							{//check if a potential increase in volume and weight is allowed on each tour
									
								double old_d = w->con[x->index]->determin + x->con[y->index]->determin;
								double new_d = w->con[b->index]->determin + b->con[y->index]->determin;
								double old_e = a->con[b->index]->determin + b->con[c->index]->determin;
								double new_e = a->con[x->index]->determin + x->con[c->index]->determin;
								double lb = (old_d - new_d) + (old_e - new_e); // positive is good
								if (lb > 1e-9)
								{
									//local evaluation on path d
									
									//the break can come too early due to potential traveltime decrease 
									bool reqbreakd = false;
									int breakindexd = tourd->breakindex;
									if (tourd->breakindex >= i + 1)//update break if it is currently position at x or after
									{
										reqbreakd = true;
									}
									double departuretime = ins->t[d].EDT + tourd->deptime[i];
									// w to b
									double arrivaltime = ins->arrival_time(w->con[b->index], departuretime);
									if ((reqbreakd) && (max(b->LTW[d] - ins->breakdur, arrivaltime) >= ins->breakstart + 1e-9))
									{
										if (max(b->LTW[d] - ins->breakdur, arrivaltime) > ins->breakend + 1e-9)
										{
											continue;
										}
										reqbreakd = false;
										breakindexd = i + 1;
										arrivaltime += ins->breakdur;

									}
									if (arrivaltime < b->LTW[d])
									{
										arrivaltime = b->LTW[d];
									}
									if (arrivaltime > b->UTW[d] + 1e-9)
									{
										continue;
									}
									arrivaltime += b->serv;
									departuretime = arrivaltime;
									//traveltime from b to y (can be on the end depot)
									arrivaltime = ins->arrival_time(b->con[y->index], departuretime);
									if ((reqbreakd) && ((max(y->LTW[d] - ins->breakdur, arrivaltime) >= ins->breakstart + 1e-9) || (y->index == ins->maxvertices - 1)))
									{
										reqbreakd = false;
										arrivaltime += ins->breakdur;
										breakindexd = i + 2;
									}
									if (arrivaltime < y->LTW[d])
									{
										arrivaltime = y->LTW[d];
									}
									if (arrivaltime > y->UTW[d] + 1e-9)
									{
										continue;
									}
									double arrivaltimey = arrivaltime + y->serv;
									double diffd = (ins->t[d].EDT + tourd->deptime[i + 2]) - arrivaltimey;
									//local evaluation on path e
									//the break can come too early due to potential traveltime decrease 
									bool reqbreake = false;
									int breakindexe = toure->breakindex;
									if (toure->breakindex >= j + 1)
									{
										reqbreake = true;
									}
									departuretime = ins->t[e].EDT + toure->deptime[j];
									//traveltime a to x
									arrivaltime = ins->arrival_time(a->con[x->index], departuretime);
									if ((reqbreake) && (max(x->LTW[e] - ins->breakdur, arrivaltime) >= ins->breakstart + 1e-9))
									{
										if (max(x->LTW[e] - ins->breakdur, arrivaltime) > ins->breakend + 1e-9)
										{
											continue;
										}
										reqbreake = false;
										breakindexe = j + 1;
										arrivaltime += ins->breakdur;
									}
									if (arrivaltime < x->LTW[e])
									{
										arrivaltime = x->LTW[e];
									}
									if (arrivaltime > x->UTW[e] + 1e-9)
									{
										continue;
									}
									arrivaltime += x->serv;
									departuretime = arrivaltime;
									//traveltime from x to c (can be the end depot)
									arrivaltime = ins->arrival_time(x->con[c->index], departuretime);
									if ((reqbreake) && ((max(c->LTW[e] - ins->breakdur, arrivaltime) >= ins->breakstart + 1e-9) || (c->index == ins->maxvertices - 1)))
									{
										reqbreake = false;
										arrivaltime += ins->breakdur;
										breakindexe = j + 2;
									}
									if (arrivaltime < c->LTW[e])
									{
										arrivaltime = c->LTW[e];
									}
									if (arrivaltime > c->UTW[e] + 1e-9)
									{
										continue;
									}
									double arrivaltimec = arrivaltime + c->serv;
									double diffe = (ins->t[e].EDT + toure->deptime[j + 2]) - arrivaltimec;
									double localdecreasetotal = diffd + diffe;
									//local travel time gain
									if (localdecreasetotal > bestdecrease + 1e-9)
									{
										//check enddepot time on path d
										bool feasible_d = true;
										double departuretime_d = arrivaltimey;
										double arrivaltime_d = arrivaltimey;
										for (int v = i + 3; v < int(tourd->seq.size()); ++v)
										{
											Ins::Vertex* o = tourd->seq[v - 1];
											Ins::Vertex* p = tourd->seq[v];
											arrivaltime_d = ins->arrival_time(o->con[p->index], departuretime_d);
											if ((reqbreakd) && ((max(p->LTW[d] - ins->breakdur, arrivaltime_d) >= ins->breakstart + 1e-9) || (p->index == ins->maxvertices - 1)))
											{
												if (arrivaltime_d > ins->breakend + 1e-9)
												{
													feasible_d = false;
													break;
												}
												breakindexd = v;
												reqbreakd = false;
												arrivaltime_d += ins->breakdur;
											}
											if (arrivaltime_d < p->LTW[d])
											{
												arrivaltime_d = p->LTW[d];
											}
											if (arrivaltime_d > p->UTW[d] + 1e-9)
											{
												feasible_d = false;
												break;
											}
											arrivaltime_d += p->serv;
											departuretime_d = arrivaltime_d;
										}
										
										//check enddepot time on path e
										bool feasible_e = true;
										double departuretime_e = arrivaltimec;
										double arrivaltime_e = arrivaltimec;
										for (int v = j + 3; v < int(toure->seq.size()); ++v)
										{
											Ins::Vertex* o = toure->seq[v - 1];
											Ins::Vertex* p = toure->seq[v];
											arrivaltime_e = ins->arrival_time(o->con[p->index], departuretime_e);
											if ((reqbreake) && ((max(p->LTW[e] - ins->breakdur, arrivaltime_e) >= ins->breakstart + 1e-9) || (p->index == ins->maxvertices - 1)))
											{
												if (arrivaltime_e > ins->breakend + 1e-9)
												{
													feasible_e = false;
													break;
												}
												breakindexe = v;
												reqbreake = false;
												arrivaltime_e += ins->breakdur;
											}
											if (arrivaltime_e < p->LTW[e])
											{
												arrivaltime_e = p->LTW[e];
											}
											if (arrivaltime_e > p->UTW[e] + 1e-9)
											{
												feasible_e = false;
												break;
											}
											arrivaltime_e += p->serv;
											departuretime_e = arrivaltime_e;
										}
										if ((!feasible_d) || (!feasible_e))
										{
											continue;   // skip the rest of this candidate (go to next j)
										}
										double globaldecreasetotal =(tourd->deptime.back() - (arrivaltime_d - ins->t[d].EDT)) +(toure->deptime.back() - (arrivaltime_e - ins->t[e].EDT));
										//global improvement check
										if (globaldecreasetotal > bestdecrease + 1e-9)
										{
											bestdecrease = globaldecreasetotal;
											beste = toure;
											bestj = j + 1;
											bestd = tourd;
											besti = i + 1;
											bestbreakindexd = breakindexd;
											bestbreakindexe = breakindexe;
											improvement = true;
											if (mode == 0)
											{
												goto swap2;
											}
										}
									}//end local check
								}//end tid check
							}//end if still capacity free
						}//end if neighbour
					}//end for j
				}//end for path e
			}//end for master sol i
		}// end for masterpath d
		if (improvement)
		{
		swap2:
			Sol remember = sol;
			Ins::Vertex *x = bestd->seq[besti];
			Ins::Vertex *b = beste->seq[bestj];
			//sol.replace_vertex(*bestd, b, besti);
			//sol.replace_vertex(*beste, x, bestj);
			sol.replace_vertex(*bestd,b,besti,bestbreakindexd);
			sol.replace_vertex(*beste,x,bestj,bestbreakindexe);
			sol.available[bestd->seq[besti]->index] = false;
			sol.available[beste->seq[bestj]->index] = false;
			double actualdecrease = 0.0;
			for (int t = 0; t < (int) sol.tours.size(); ++t)
			{
				actualdecrease += remember.tours[t].deptime.back()-sol.tours[t].deptime.back();
			}
			if (!sol.check())
			{
				cout<< term::fg(term::Color::red) << "sol error in swap2" << endl;
			}
			if (fabs(bestdecrease - actualdecrease) > 0.01)
			{
				cout<< term::fg(term::Color::red) << "time decrease error swap2: "<<bestdecrease-actualdecrease << endl;
			}
			succes = true;
		}//end if improvement
	}//end while improvement
	return succes;
}//end swap2
