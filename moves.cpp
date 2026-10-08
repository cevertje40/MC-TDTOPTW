#include "moves.h"

using namespace std;


void Moves::parallel_construct(Sol& sol)
{
	// 0) Per-tour depot-reachable lists, sorted by score
	std::vector<std::vector<Ins::Vertex*>> dep_reach(ins->maxtours);
	for (int d = 0; d < ins->maxtours; ++d)
	{
		for (int i = 1; i < ins->maxvertices - 1; ++i)
		{
			const double lb = ins->v[0].con[i]->determin + ins->v[i].serv + ins->v[i].con[ins->maxvertices - 1]->determin;
			if (lb <= ins->t[d].T_max + 1e-9)
				dep_reach[d].push_back(&ins->v[i]);
		}
		std::sort(dep_reach[d].begin(), dep_reach[d].end(), [](Ins::Vertex* a, Ins::Vertex* b) { return a->score > b->score; });
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

			double at = ins->arrival_time(last->con[cand->index], currenttime);
			int need_break = 0;
			if (tour.breakindex == -1 && at >= ins->breakstart - 1e-9)
			{
				// Arrival is in the break window: take the break at this customer.
				if (at > ins->breakend + 1e-9) return false;

				at += ins->breakdur;
				need_break = 1;
			}

			// Apply the customer's time window after any break.
			if (at < cand->LTW[d]) at = cand->LTW[d];
			if (at > cand->UTW[d] + 1e-9) return false;

			at += cand->serv;

			double enddep =ins->arrival_time(cand->con[ins->maxvertices - 1], at);

			if (tour.breakindex == -1 && need_break == 0)
			{
				// The break would be taken upon arrival at the end depot.
				if (enddep > ins->breakend + 1e-9) return false;

				enddep += ins->breakdur;
			}

			if (enddep > ins->t[d].LAT + 1e-9) return false;

			new_deptime = at - ins->t[d].EDT; // includes break if taken
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

template<RatioFn RATIO>boost::heap::priority_queue<One_one_rep_nb>Moves::one_one_replace_gen_nb_kernel(Sol& sol, TabuVector& tabulist, int globalbest)
{
	boost::heap::priority_queue<One_one_rep_nb> adm_nb;

	const int max_threads =
	#ifdef _OPENMP
		omp_get_max_threads();
	#else
		1;
	#endif

	std::vector<boost::heap::priority_queue<One_one_rep_nb>> tls_queues(max_threads);

	#pragma omp parallel
	{
		const int tid =
		#ifdef _OPENMP
			omp_get_thread_num();
		#else
			0;
		#endif

		auto& local_q = tls_queues[tid];

		#pragma omp for nowait
		for (int d = 0; d < ins->maxtours; ++d)
		{
			double local_bestkey = -DBL_MAX;
			// ===== INSERT PART =====
			Sol::Tour& tour = sol.tours[d];
			int endj = (int)tour.seq.size();
			std::vector<double> wait_suffix = tour.compute_wait_suffix();
			const auto& Td = ins->t[d];
			const double EDT = Td.EDT;
			const double Wmax = Td.W_max, Vmax = Td.V_max, Tmax = Td.T_max;
			const double breakdur = ins->breakdur;
			for (int j = 0; j < endj - 1; ++j)
			{
				Ins::Vertex* x = tour.seq[j];
				Ins::Vertex* z = tour.seq[j + 1];
				int breakz = tour.action[j + 1];
				int nb_size = (int)x->nb[d].size();

				double currenttime = tour.deptime[j] + EDT;

				for (int i = 0; i < nb_size - 1; ++i) // skip depot
				{
					Ins::Vertex* y = x->nb[d][i];

					if (sol.score + y->score <= globalbest + 1e-9)
					{
						if (tabulist.isTabu(y->index, d)) continue;
					}

					if (!(sol.available[y->index] && x->nbi[d][y->index] && y->nbi[d][z->index])) continue;

						if (tour.weight + y->weight > Wmax + 1e-9) continue;
						if (tour.volume + y->volume > Vmax + 1e-9) continue;

						double t0 = tour.deptime[j] + EDT;
						double base_up = ins->arrival_time(x->con[z->index], t0) - t0;
						double ins_lb = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
						double dt_lb = std::max(0.0, ins_lb - base_up);

					double key_ub = RATIO(dt_lb, y->score, y->weight, y->volume,Tmax, Wmax, Vmax);

					if ((key_ub <= local_bestkey + 1e-9) || (dt_lb > wait_suffix[j + 1] + tour.max_shift[j + 1] + 1e-9)) continue;

					double at = ins->arrival_time(x->con[y->index], currenttime);

					if (at < y->LTW[d]) at = y->LTW[d];

					if (at > y->UTW[d] + 1e-9) continue;

					at += y->serv;

					at = ins->arrival_time(y->con[z->index], at);

					if (breakz)
					{
						const bool endp = (z->index == ins->maxvertices - 1);
						// - Non-depot: if arrive early, wait to breakstart, then add B
						// - Depot: allowed to start upon arrival even if < breakstart
						if (!endp && at < ins->breakstart) at = ins->breakstart;
						// Must still start the break not after breakend
						if (at > ins->breakend + 1e-9) continue;  // no legal break start here
						at += breakdur;  // finish break at z
					}

					if (at < z->LTW[d]) at = z->LTW[d];
					at += z->serv;
					
					double shift = (at - EDT) - tour.deptime[j + 1];
					if (shift > tour.max_shift[j + 1] + 1e-9) continue;

					double key = RATIO(shift, y->score, y->weight, y->volume, Tmax, Wmax, Vmax);

					if (key > local_bestkey + 1e-9)
					{
						local_bestkey = key;
						local_q.push(One_one_rep_nb(d, -1, j, y, y->score, key, { {},{y} }));
					}
				}
			}

			// ===== REPLACE PART =====
			int endh = (int)sol.tours[d].seq.size();
			for (int h = 1; h < endh - 1; ++h)
			{
				Tour tourrem = sol.tours[d];
				Ins::Vertex* r = tourrem.seq[h];
				tourrem.remove_vertex(h);
				std::vector<double> wait_suffix_r = tourrem.compute_wait_suffix();

				int endj2 = (int)tourrem.seq.size();
				for (int j = 0; j < endj2 - 1; ++j)
				{
					Ins::Vertex* x = tourrem.seq[j];
					Ins::Vertex* z = tourrem.seq[j + 1];
					int breakz = tourrem.action[j + 1];
					int nb_size = (int)x->nb[d].size();

					double currenttime = tourrem.deptime[j] + EDT;

					for (int i = 0; i < nb_size - 1; ++i)
					{
						Ins::Vertex* y = x->nb[d][i];

						if (sol.score + (y->score - r->score) <= globalbest + 1e-9)
						{
							if ((tabulist.isTabu(r->index, d)) || (tabulist.isTabu(y->index, d))) continue;
						}
						
						if (!(sol.available[y->index] && x->nbi[d][y->index] && y->nbi[d][z->index])) continue;

						if (tourrem.weight + y->weight > Wmax + 1e-9) continue;
						if (tourrem.volume + y->volume > Vmax + 1e-9) continue;

						double t0      = tourrem.deptime[j] + EDT;
						double base_up = ins->arrival_time(x->con[z->index], t0) - t0;
						double ins_lb  = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
						double dt_lb   = std::max(0.0, ins_lb - base_up);

						double key_ub  = RATIO(dt_lb,y->score - r->score,y->weight - r->weight,y->volume - r->volume,Tmax, Wmax, Vmax);

						if ((key_ub <= local_bestkey + 1e-9) ||
							(dt_lb   >  wait_suffix_r[j + 1] + tourrem.max_shift[j + 1] + 1e-9))
							continue;

						double at = ins->arrival_time(x->con[y->index], currenttime);
						if (at < y->LTW[d]) at = y->LTW[d];
						if (at > y->UTW[d] + 1e-9) continue;
						at += y->serv;
						at  = ins->arrival_time(y->con[z->index], at);
						if (breakz)
						{
							const bool endp = (z->index == ins->maxvertices - 1);
							// - Non-depot: if arrive early, wait to breakstart, then add B
							// - Depot: allowed to start upon arrival even if < breakstart
							if (!endp && at < ins->breakstart) at = ins->breakstart;
							// Must still start the break not after breakend
							if (at > ins->breakend + 1e-9) continue;  // no legal break start here
							at += breakdur;  // finish break at z
						}
						if (at < z->LTW[d]) at = z->LTW[d];
						at += z->serv;
						
						double shift = (at - EDT) - tourrem.deptime[j + 1];
						if (shift > tourrem.max_shift[j + 1] + 1e-9) continue;

						double key = RATIO(shift,y->score - r->score,y->weight - r->weight,y->volume - r->volume,Tmax, Wmax, Vmax);

						if (key > local_bestkey + 1e-9) 
						{
							local_bestkey = key;
							local_q.push(One_one_rep_nb(d, h, j, y,y->score - r->score, key, { {r},{y} }));
						}
					}
				}
			}
		} // for d
	} // parallel

	// Serial merge (outside hot path)
	for (auto& q : tls_queues) 
	{
		while (!q.empty()) 
		{
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

template<RatioFn RATIO>boost::heap::priority_queue<Two_one_rep_nb>Moves::two_one_replace_gen_nb_kernel(Sol& sol, TabuVector& tabulist, int globalbest)
{
	boost::heap::priority_queue<Two_one_rep_nb> adm_nb;

	const int max_threads =
	#ifdef _OPENMP
		omp_get_max_threads();
	#else
		1;
	#endif

	std::vector<boost::heap::priority_queue<Two_one_rep_nb>> tls_queues(max_threads);

	#pragma omp parallel
	{
		const int tid =
		#ifdef _OPENMP
			omp_get_thread_num();
		#else
			0;
		#endif

		auto& local_q = tls_queues[tid];

		// ----------------------- PARALLEL BODY -----------------------
		#pragma omp for nowait
		for (int d = 0; d < ins->maxtours; ++d)
		{
			// IMPORTANT: reset per tour d (matches your original behavior)
			double local_bestkey = -DBL_MAX;
			const auto& Td = ins->t[d];
			const double EDT = Td.EDT;
			const double Wmax = Td.W_max, Vmax = Td.V_max, Tmax = Td.T_max;
			const double breakdur = ins->breakdur;
			int end = (int)sol.tours[d].seq.size();
			for (int g = 1; g < end - 1; ++g)            // first removed regular vertex
			{
				for (int h = g + 1; h < end - 1; ++h)    // second removed regular vertex
				{
					// Build tour with 2 removals
					Sol::Tour tourrem = sol.tours[d];
					Ins::Vertex* r = tourrem.seq[g];
					Ins::Vertex* s = tourrem.seq[h];
					tourrem.remove_vertices(g, h);

					std::vector<double> wait_suffix = tourrem.compute_wait_suffix();

					const int endj = (int)tourrem.seq.size();
					for (int j = 0; j < endj - 1; ++j) // insertion boundary in tourrem
					{
						Ins::Vertex* x = tourrem.seq[j];
						Ins::Vertex* z = tourrem.seq[j + 1];
						const int breakz = tourrem.action[j + 1];

						// Hoist base values that don't depend on y:
						const double t0 = tourrem.deptime[j] + EDT;
						const double base_up = ins->arrival_time(x->con[z->index], t0) - t0;

						const int nb_size = (int)x->nb[d].size();
						const double currenttime = t0; 

						for (int i = 0; i < nb_size - 1; ++i) // skip last neighbor = end depot
						{
							Ins::Vertex* y = x->nb[d][i];

							// Tabu with aspiration against globalbest
							if ((sol.score + (y->score - (r->score + s->score)) <= globalbest + 1e-9))
							{
								if ((tabulist.isTabu(r->index, d)) ||(tabulist.isTabu(s->index, d)) ||	(tabulist.isTabu(y->index, d)))	continue;
								
							}
							
							// Availability & adjacency
							if (!(sol.available[y->index] && x->nbi[d][y->index] && y->nbi[d][z->index]))
								continue;

							// Capacity
							if (tourrem.weight + y->weight > Wmax + 1e-9) continue;
							if (tourrem.volume + y->volume > Vmax + 1e-9) continue;

							// Time LB/UB gate (cheap)
							const double ins_lb = x->con[y->index]->determin + y->serv + y->con[z->index]->determin;
							const double dt_lb = std::max(0.0, ins_lb - base_up);

							const double key_ub = RATIO(dt_lb,y->score - (r->score + s->score),y->weight - (r->weight + s->weight),y->volume - (r->volume + s->volume),ins->t[d].T_max, ins->t[d].W_max, ins->t[d].V_max);

							if ((key_ub <= local_bestkey + 1e-9) ||	(dt_lb > wait_suffix[j + 1] + tourrem.max_shift[j + 1] + 1e-9))
							{
								continue; // prune
							}

							// Exact time-feasibility check
							double at = ins->arrival_time(x->con[y->index], currenttime);
							if (at < y->LTW[d]) at = y->LTW[d];
							if (at > y->UTW[d] + 1e-9) continue;
							at += y->serv;

							at = ins->arrival_time(y->con[z->index], at);
							if (breakz) 
							{
								const bool endp = (z->index == ins->maxvertices - 1);

								// - Non-depot: if arrive early, wait to breakstart, then add B
								// - Depot: allowed to start upon arrival even if < breakstart
								if (!endp && at < ins->breakstart) at = ins->breakstart;
								// Must still start the break not after breakend
								if (at > ins->breakend + 1e-9) continue;  // no legal break start here
								at += breakdur;  // finish break at z
							}

							// LTW at z is enforced AFTER the break (or immediately if no break at z)
							if (at < z->LTW[d]) at = z->LTW[d];
							at += z->serv;

							const double shift = (at - EDT) - tourrem.deptime[j + 1];
							if (shift > tourrem.max_shift[j + 1] + 1e-9) continue;

							// Objective key
							const double key = RATIO(shift,y->score - (r->score + s->score),y->weight - (r->weight + s->weight),y->volume - (r->volume + s->volume),Tmax, Wmax, Vmax);

							if (key > local_bestkey + 1e-9) 
							{
								local_bestkey = key;
								// Push to this thread's queue (no locking)
								local_q.push(Two_one_rep_nb(d, g, h, j, y,y->score - (r->score + s->score),key, { {r, s}, {y}}));
							}
						} // nb
					} // j
				} // h
			} // g
		} // for d
		// --------------------- END PARALLEL BODY ---------------------
	} // omp parallel

	// Serial merge 
	for (auto& q : tls_queues) 
	{
		while (!q.empty()) 
		{
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

template<RatioFn RATIO>boost::heap::priority_queue<One_two_rep_nb>Moves::one_two_replace_gen_nb_kernel(Sol& sol, TabuVector& tabulist, int globalbest)
{
	boost::heap::priority_queue<One_two_rep_nb> adm_nb;

	const int max_threads =
	#ifdef _OPENMP
		omp_get_max_threads();
	#else
		1;
	#endif

	std::vector<boost::heap::priority_queue<One_two_rep_nb>> tls_queues(max_threads);

	#pragma omp parallel
	{
		const int tid =
		#ifdef _OPENMP
			omp_get_thread_num();
		#else
			0;
		#endif

		auto& local_q = tls_queues[tid];

		// ----------------------- PARALLEL BODY -----------------------
		#pragma omp for schedule(guided)
		for (int d = 0; d < ins->maxtours; ++d)
		{
			const auto& Td = ins->t[d];
			const double EDT = Td.EDT;
			const double Wmax = Td.W_max, Vmax = Td.V_max, Tmax = Td.T_max;
			const double breakdur = ins->breakdur;
			// reset per tour (matches your original behavior)
			double local_bestkey = -DBL_MAX;

			int endh = (int)sol.tours[d].seq.size();
			for (int h = 1; h < endh - 1; ++h) // removed vertex r
			{
				// build tour after removing r
				Sol::Tour tourrem = sol.tours[d];
				Ins::Vertex* r = tourrem.seq[h];
				tourrem.remove_vertex(h);
				const double W_free_after_rem = Wmax - (sol.tours[d].weight - r->weight);
				const double V_free_after_rem = Vmax - (sol.tours[d].volume - r->volume);
				// suffix waiting (upper bound) on tourrem
				std::vector<double> wait_suffix_rem = tourrem.compute_wait_suffix();

				int endj = (int)tourrem.seq.size();
				for (int j = 0; j < endj - 1; ++j) // position for first insertion y1
				{
					Ins::Vertex* x1 = tourrem.seq[j];
					Ins::Vertex* z1 = tourrem.seq[j + 1];
					int breakz1 = tourrem.action[j + 1];

					int nb1 = (int)x1->nb[d].size();
					for (int i = 0; i < nb1 - 1; ++i) // skip end depot
					{
						Ins::Vertex* y1 = x1->nb[d][i];

						if (!(sol.available[y1->index] && x1->nbi[d][y1->index] && y1->nbi[d][z1->index]))
							continue;

						// capacity with y1 (remember we already removed r)
						if (y1->weight > W_free_after_rem + 1e-9) continue;
						if (y1->volume > V_free_after_rem + 1e-9) continue;

						// time LB/UB gate for first insertion
						double t0_1 = tourrem.deptime[j] + EDT;
						double baseUB_1 = ins->arrival_time(x1->con[z1->index], t0_1) - t0_1;
						double insLB_1 = x1->con[y1->index]->determin + y1->serv + y1->con[z1->index]->determin;
						double dtLB_1 = std::max(0.0, insLB_1 - baseUB_1);

						// feasibility pre-gate using waiting budget on tourrem
						if (dtLB_1 > wait_suffix_rem[j + 1] + tourrem.max_shift[j + 1] + 1e-9)
							continue;

						// exact feasibility for first insertion (shift1)
						double ct = t0_1;
						double at = ins->arrival_time(x1->con[y1->index], ct);
						if (at < y1->LTW[d]) at = y1->LTW[d];
						if (at > y1->UTW[d] + 1e-9) continue;
						at += y1->serv;

						at = ins->arrival_time(y1->con[z1->index], at);
						if (breakz1) 
						{
							const bool endp = (z1->index == ins->maxvertices - 1);

							// Start-only rule:
							// - Non-depot: if arrive early, wait to breakstart, then add B
							// - Depot: allowed to start upon arrival even if < breakstart
							if (!endp && at < ins->breakstart) at = ins->breakstart;

							// Must still start the break not after breakend
							if (at > ins->breakend + 1e-9) continue;  // no legal break start here

							at += breakdur;  // finish break at z
						}

						// LTW at z is enforced AFTER the break (or immediately if no break at z)
						if (at < z1->LTW[d]) at = z1->LTW[d];
						at += z1->serv;

						double shift1 = (at - EDT) - tourrem.deptime[j + 1];
						if (shift1 > tourrem.max_shift[j + 1] + 1e-9) continue;

						// ---------------- Stage-A fast pre-gate for second insertion ----------------
						bool any_second_survives = false;
						const double W_free_after_y1 = W_free_after_rem - y1->weight;
						const double V_free_after_y1 = V_free_after_rem - y1->volume;
						// boundaries p != j in tourrem (map to k != j after insertion)
						int p_end = (int)tourrem.seq.size() - 1;
						for (int p = 0; p < p_end && !any_second_survives; ++p) {
							if (p == j) continue;

							Ins::Vertex* x2 = tourrem.seq[p];
							Ins::Vertex* z2 = tourrem.seq[p + 1];

							// delay reaching boundary p+1 after first insertion (UB)
							double Wprefix = 0.0;
							if (p + 1 > j + 1) {
								Wprefix = wait_suffix_rem[j + 1] - wait_suffix_rem[p + 1];
								if (Wprefix < 0.0) Wprefix = 0.0;
							}
							double delta_at_p = shift1 - Wprefix;
							if (delta_at_p < 0.0) delta_at_p = 0.0;

							double t_low = tourrem.deptime[p] + EDT;
							double t_high = t_low + delta_at_p;

							// base-arc travel-time UB over [t_low, t_high]
							Ins::Connec* cxz = x2->con[z2->index];
							double baseUB = base_travel_UB_over_interval(ins, cxz, t_low, t_high);

							double Wup = wait_suffix_rem[p + 1];
							double Mup = tourrem.max_shift[p + 1];

							// neighbours from x2 (skip end depot)
							std::vector<Ins::Vertex*>& neigh2 = x2->nb[d];
							int nb2 = (int)neigh2.size();
							for (int ii = 0; ii < nb2 - 1; ++ii) {
								Ins::Vertex* y2 = neigh2[ii];

								// optional tabu
								if ((sol.score + (y1->score + y2->score - r->score) <= globalbest + 1e-9))
								{
									if ((tabulist.isTabu(r->index, d)) || (tabulist.isTabu(y1->index, d)) || (tabulist.isTabu(y2->index, d))) continue;
								}
								 
								// capacity with y1 + y2
								if (y2->weight > W_free_after_y1 + 1e-9) continue;
								if (y2->volume > V_free_after_y1 + 1e-9) continue;

								// LB on inserted path at boundary p
								double insLB = x2->con[y2->index]->determin + y2->serv + y2->con[z2->index]->determin;
								double dtLB = insLB - baseUB; if (dtLB < 0.0) dtLB = 0.0;

								if (dtLB <= Wup + Mup + 1e-9) { any_second_survives = true; break; }
							}
						}

						// special-case boundary k == j (between x1 and y1 in tourremins)
						if (!any_second_survives) {
							Ins::Connec* c_x1y1 = x1->con[y1->index];
							double baseUB_j = ins->arrival_time(c_x1y1, t0_1) - t0_1;

							double Wup_j = wait_suffix_rem[j + 1];
							double Mup_j = tourrem.max_shift[j + 1];

							std::vector<Ins::Vertex*>& neigh1 = x1->nb[d];
							int nbx = (int)neigh1.size();
							for (int ii = 0; ii < nbx - 1; ++ii) {
								Ins::Vertex* y2 = neigh1[ii];
								if (y2 == y1) continue;

								if (sol.score + (y1->score + y2->score - r->score) <= globalbest + 1e-9)
								{
									if ((tabulist.isTabu(r->index, d)) || (tabulist.isTabu(y1->index, d)) || (tabulist.isTabu(y2->index, d))) continue;
								}

								if (y2->weight > W_free_after_y1 + 1e-9) continue;
								if (y2->volume > V_free_after_y1 + 1e-9) continue;

								double insLB = x1->con[y2->index]->determin + y2->serv + y2->con[y1->index]->determin;
								double dtLB = insLB - baseUB_j; if (dtLB < 0.0) dtLB = 0.0;

								if (dtLB <= Wup_j + Mup_j + 1e-9) { any_second_survives = true; break; }
							}
						}

						if (!any_second_survives) continue; // skip (j,y1) if no second insertion can survive

						// ---------------- Do the actual first insertion; then test all k,y2 ----------------
						Sol::Tour tourremins = tourrem;
						tourremins.insert_vertex(y1, j);

						std::vector<double> wait_suffix = tourremins.compute_wait_suffix();
						int endk = (int)tourremins.seq.size();

						for (int k = 0; k < endk - 1; ++k)
						{
							Ins::Vertex* x2 = tourremins.seq[k];
							Ins::Vertex* z2 = tourremins.seq[k + 1];
							int breakz2 = tourremins.action[k + 1];

							const double t0 = tourremins.deptime[k] + EDT;
							const double base_up = ins->arrival_time(x2->con[z2->index], t0) - t0;

							int nb2 = (int)x2->nb[d].size();
							for (int ii = 0; ii < nb2 - 1; ++ii) // skip end depot
							{
								Ins::Vertex* y2 = x2->nb[d][ii];

								if (sol.score + (y1->score + y2->score - r->score) <= globalbest + 1e-9)
								{
									if ((tabulist.isTabu(r->index, d)) || (tabulist.isTabu(y1->index, d)) || (tabulist.isTabu(y2->index, d))) continue;
								}

								if (!(sol.available[y2->index] && x2->nbi[d][y2->index] && y2->nbi[d][z2->index]) || (y2 == y1))
									continue;

								if (y2->weight > W_free_after_y1 + 1e-9) continue;
								if (y2->volume > V_free_after_y1 + 1e-9) continue;

								// cheap key_ub + shift gate
								double ins_lb = x2->con[y2->index]->determin + y2->serv + y2->con[z2->index]->determin;
								double dt_lb = std::max(0.0, ins_lb - base_up);
								double key_ub = RATIO(shift1 + dt_lb,(y1->score + y2->score) - r->score,(y1->weight + y2->weight) - r->weight,(y1->volume + y2->volume) - r->volume,Tmax, Wmax, Vmax);

								if (key_ub <= local_bestkey + 1e-9) continue;
								if (dt_lb > wait_suffix[k + 1] + tourremins.max_shift[k + 1] + 1e-9) continue;

								// exact second insertion (shift2)
								double ct2 = t0;
								double at2 = ins->arrival_time(x2->con[y2->index], ct2);
								if (at2 < y2->LTW[d]) at2 = y2->LTW[d];
								if (at2 > y2->UTW[d] + 1e-9) continue;
								at2 += y2->serv;

								at2 = ins->arrival_time(y2->con[z2->index], at2);
								if (breakz2) 
								{
									const bool endp = (z2->index == ins->maxvertices - 1);

									// - Non-depot: if arrive early, wait to breakstart, then add B
									// - Depot: allowed to start upon arrival even if < breakstart
									if (!endp && at2 < ins->breakstart) at2 = ins->breakstart;

									// Must still start the break not after breakend
									if (at2 > ins->breakend + 1e-9) continue;  // no legal break start here

									at2 += breakdur;  // finish break at z
								}

								// LTW at z is enforced AFTER the break (or immediately if no break at z)
								if (at2 < z2->LTW[d]) at2 = z2->LTW[d];
								at2 += z2->serv;

								double shift2 = (at2 - EDT) - tourremins.deptime[k + 1];
								if (shift2 > tourremins.max_shift[k + 1] + 1e-9) continue;

								// final key
								double key = RATIO(shift1 + shift2,(y1->score + y2->score) - r->score,(y1->weight + y2->weight) - r->weight,(y1->volume + y2->volume) - r->volume,Tmax, Wmax, Vmax);

								if (key > local_bestkey + 1e-9) 
								{
									local_bestkey = key;
									local_q.push(One_two_rep_nb(d, h, j, k, y1, y2,(y1->score + y2->score) - r->score,key,{ {r}, {y1, y2} }
									));
								}
							} // y2
						} // k
					} // y1
				} // j
			} // h
		} // for d
		// --------------------- END PARALLEL BODY ---------------------
	} // omp parallel

	// Serial merge (outside hot path)
	for (auto& q : tls_queues) {
		while (!q.empty()) {
			adm_nb.push(q.top()); // copy
			q.pop();
		}
	}

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
								//do not account for break repositioning so don't use linker on action
								if (arrivaltime + tour.action[l + 1] * (ins->breakdur) < second->LTW[d])
								{
									arrivaltime = second->LTW[d] - (tour.action[l + 1] * ins->breakdur);
								}
								if (arrivaltime > second->UTW[d] + 1e-9)
								{
									currenttime = DBL_MAX;
									break;
								}
								arrivaltime += second->serv + (tour.action[l + 1] * ins->breakdur);
								currenttime = arrivaltime;
							}// end l
							newtraveltime = currenttime - departuretime;
							oldtraveltime = (ins->t[d].EDT + tour.deptime[j + 1]) - departuretime;
							delta_tt = oldtraveltime - newtraveltime;
							if (delta_tt > bestdelta+ 1e-9)//local evaluation
							{
								bool reqbreak = true;
								currenttime = ins->t[d].EDT;
								for (int l = 0; l < end - 1; ++l)//global evaluation
								{
									Ins::Vertex* first = tour.seq[map_idx(l)];
									Ins::Vertex* second = tour.seq[map_idx(l + 1)];//can be the end depot
									double arrivaltime = ins->arrival_time(first->con[second->index], currenttime);
									//account for break
									if ((reqbreak) && ((max(second->LTW[d] - ins->breakdur,arrivaltime) >= ins->breakstart + 1e-9) || (second->index == ins->maxvertices - 1)))
									{
										arrivaltime += ins->breakdur;
										reqbreak = false;
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
								if (reqbreak == true)//als je aankomt bij het einddepot en nog steeds geen break genomen hebt
								{
									currenttime += ins->breakdur;//breaktime bijtellen bij aankomst tijd bij einddepot
								}
								newtraveltime = currenttime - departuretime;
								oldtraveltime = ((ins->t[d].EDT + tour.deptime[end - 1]) - departuretime);
								delta_tt = oldtraveltime-newtraveltime;
								if (delta_tt >= bestdelta + 1e-9)//global evaluation
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
						}//end tid filter
					}//end neighbourhood check
				}//end j
			}// end i
			if (improvement)
			{
				swap:
				//Sol remember = sol;
				tour.swap_vertices(besti, bestj);
				/*
				double actualdecrease = remember.tours[d].deptime.back() - tour.deptime.back();
				if (!sol.check())
				{
					cout << "error in swap" << endl;
				}
				if (fabs(bestdelta - actualdecrease) > 0.01)
				{
					cout << "error swap" << endl;
				}
				*/
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
							//int remember = sol.traveltime[d].back();
							if (currenttime == DBL_MAX) continue;
							delta_tt = (ins->t[d].EDT + tour.deptime[j + 1]) - currenttime;
							if (delta_tt > +1e-9)//local evaluation
							{
								Tour temptour = tour;
								temptour.opt_vertices(i, j);
								delta_tt = tour.deptime.back() - temptour.deptime.back();
								if (delta_tt > bestdelta + 1e-9)//global evaluation
								{
									improvement = true;
									bestdelta = delta_tt;
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
				//Sol remember = sol;
				tour.opt_vertices(besti,bestj);
				//double actualdecrease = remember.tours[d].deptime.back() - tour.deptime.back();
				/*
				if (!sol.check())
				{
					cout << "error in two-opt" << endl;
				}
				if (fabs(bestdelta - actualdecrease) > 0.01)
				{
					cout << "error in two-opt" << endl;
				}
				*/
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
				double ttwxy = (tourd->deptime[i + 2] - (y->serv + tourd->action[i + 2] * ins->breakdur + x->serv + tourd->action[i + 1] * ins->breakdur)) - tourd->deptime[i];
				for (int e = 0; e < ins->maxtours; ++e)
				{
					if (d != e)
					{
						Sol::Tour* toure = &sol.tours[e];
						for (int j = 0; j < int(toure->seq.size()) - 2; ++j)
						{
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
										double ttab = (toure->deptime[j + 1] - (b->serv + toure->action[j + 1] * ins->breakdur)) - toure->deptime[j];
										//calculate axb
										double departuretime = ins->t[e].EDT + toure->deptime[j];
										//traveltime a to x
										double arrivaltime = ins->arrival_time(a->con[x->index], departuretime);
										if (arrivaltime < x->LTW[e])
										{
											arrivaltime = x->LTW[e];
										}
										if (arrivaltime > x->UTW[e] + 1e-9)
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
										if (arrivaltime > b->UTW[e] + 1e-9)
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
										if ((reqbreak) && ((max(y->LTW[d] - ins->breakdur, arrivaltime) >= ins->breakstart + 1e-9) || (y->index == ins->maxvertices - 1)))
										{//break na y wordt bij global evaluation in rekening gebracht
											if (arrivaltime > ins->breakend + 1e-9)//pushing to break to right after removal does not work
											{
												continue;
											}
											reqbreak = false;
											arrivaltime += ins->breakdur;
										}
										if (arrivaltime < y->LTW[d])
										{
											arrivaltime = y->LTW[d];
										}
										if (arrivaltime > y->UTW[d] + 1e-9)
										{
											continue;
										}
										double ttwy = arrivaltime - departuretime;
										//double arrivaly = arrivaltime + y->serv;
										double shift = (arrivalb - ins->t[e].EDT) - toure->deptime[j + 1];
										double localdecreasetotal = ttwxy + ttab - (ttwy + ttaxb);
										if ((shift <= toure->max_shift[j + 1] + 1e-9) && (localdecreasetotal > bestdecrease + 1e-9))//local improvement check
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
			//Sol remember = sol;
			Ins::Vertex* candidate = bestd->seq[besti];
			sol.remove_vertex(*bestd, besti);//remove vertex from path d
			sol.insert_vertex(*beste, candidate, bestj);//insert vertex on path e
			/*
			double actualdecrease = 0.0;
			for (int t = 0; t < (int)sol.tours.size(); ++t)
			{
				actualdecrease += remember.tours[t].deptime.back() - sol.tours[t].deptime.back();
			}
			if (fabs(bestdecrease - actualdecrease) > 0.01)
			{
				cout << "error move" << endl;
			}
			*/
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
			//if (!sol.check())
			//{
				//cout << "error in move_nb" << endl;
			//}
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
			//Sol remember = sol;
			Ins::Vertex *x = bestd->seq[besti];
			Ins::Vertex *b = beste->seq[bestj];
			sol.replace_vertex(*bestd, b, besti);
			sol.replace_vertex(*beste, x, bestj);
			//sol.replace_vertex(*bestd,b,besti,bestbreakindexd);
			//sol.replace_vertex(*beste,x,bestj,bestbreakindexe);
			sol.available[bestd->seq[besti]->index] = false;
			sol.available[beste->seq[bestj]->index] = false;
			//double actualdecrease = 0.0;
			//for (int t = 0; t < (int) sol.tours.size(); ++t)
			//{
				//actualdecrease += remember.tours[t].deptime.back()-sol.tours[t].deptime.back();
			//}
			//if (!sol.check())
			//{
				//cout << "error in swap2" << endl;
			//}
			//if (fabs(bestdecrease - actualdecrease) > 0.01)
			//{
				//cout << "error swap2: "<<bestdecrease-actualdecrease << endl;
			//}
			succes = true;
		}//end if improvement
	}//end while improvement
	return succes;
}//end swap2
