#include "solution.h"

Sol::Sol(Ins& ins) :ins(&ins)
{
	solution.resize(ins.maxtours);
	traveltime.resize(ins.maxtours);
	max_shift.resize(ins.maxtours);
	scores.resize(ins.maxtours);
	weight.resize(ins.maxtours);
	volume.resize(ins.maxtours);
	action.resize(ins.maxtours);
	breakindex.resize(ins.maxtours);
	for (int t = 0; t < ins.maxtours; ++t)
	{
		solution[t].reserve(ins.maxvertices);//reserves memory
		solution[t].push_back(&ins.v[0]);//insert start depot
		traveltime[t].reserve(ins.maxvertices);//reserve memory
		traveltime[t].push_back(0);//insert first travel time
		action[t].reserve(ins.maxvertices);
		action[t].push_back(0);
		max_shift[t].reserve(ins.maxvertices);
		max_shift[t].push_back(0);
		scores[t] = 0;
		weight[t] = 0.0;
		volume[t] = 0.0;
		breakindex[t] = -1;
		tourindex.push_back(t);
	}
	shuffle(tourindex.begin(), tourindex.end(), ins.engine);
	score = 0;
	available = boost::dynamic_bitset<>(ins.maxvertices);
	available.set();//sets all bits to true
	available[ins.v[0].index] = false;
}

void Sol::check()
{
	int scorecheck = 0;
	vector<int> included(ins->maxvertices, 0);
	for (int d = 0; d < ins->maxtours; ++d)
	{
		//cout << "tour: " << d << endl << endl;
		double weightcheck = 0;
		double volumecheck = 0;
		//1. travel time check
		double currenttime = ins->t[d].EDT + traveltime[d][0];
		double length = 0;
		int end = (int)solution[d].size() - 1;
		for (int i = 0; i < end; ++i)
		{
			++included[solution[d][i]->index];
			scorecheck += solution[d][i]->score;
			weightcheck += solution[d][i]->weight;
			volumecheck += solution[d][i]->volume;
			Ins::Vertex* last = solution[d][i];
			Ins::Vertex* current = solution[d][i + 1];
			int breakcurrent = action[d][i + 1];
			double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
			double waitingtime = 0;
			if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[d])
			{
				waitingtime = current->LTW[d] - (arrivaltime + breakcurrent * (ins->breakdur));
				//cout << "waiting time for: "<<"i"<<i+1<<" , " <<sol.solution[d][i+1]->index<<" <=> "<<waitingtime << endl;
				arrivaltime = current->LTW[d] - breakcurrent * (ins->breakdur);
			}
			arrivaltime += current->serv+breakcurrent*ins->breakdur;
			//cout<<i+1<<" calc traveltime: " << arrivaltime-t[d].EDT << "stored: " << sol.traveltime[d][i+1] << endl;
			currenttime = arrivaltime;
		}//end for i
		length = currenttime - ins->t[d].EDT;
		if (length != traveltime[d].back())
			cout << red << "path: " << d << "new calculated length: " << length << " stored length: " << traveltime[d].back() << "max length" << ins->t[d].T_max << endl;
		//3. TW check
		for (int i = 0; i <= end; ++i)//utw van end depot ook checken
		{
			Ins::Vertex* current = solution[d][i];
			int breakcurrent = action[d][i];
			if (traveltime[d][i] + ins->t[d].EDT - solution[d][i]->serv < current->LTW[d])//service time zit al in traveltime
				cout << red << "path: " << d << "FAILURE!!! LTW fail for solutionnr: " << i << " /vertex index: " << current->index << endl;
			if (traveltime[d][i] + ins->t[d].EDT - solution[d][i]->serv > current->UTW[d])
			{
				cout << red << "path: " << d << "FAILURE!!! UTW fail for solutionnr: " << i << " /vertex index: " << current->index << endl;
			}
		}//end for i
		//4. weight check
		if (weightcheck != weight[d])
			cout << red << "path: " << d << "new calculated weight" << weightcheck << "stored weight: " << weight[d] << "max: " << ins->t[d].W_max << endl;
		//5. volume check
		if (volumecheck != volume[d])
			cout << yellow << "path: " << d << "new calculated volume" << volumecheck << "stored volume " << volume[d] << "max: " << ins->t[d].V_max << endl;
		//6. break timing check
		bool breakcheck = false;
		int amountbreaks = 0;
		for (int i = 0; i <= end; ++i)//break kan op enddepot zitten
		{
			if (action[d][i] == 1)
			{
				if (breakindex[d] != i)
				{
					cout << red << "path: " << d << "breakindex and sol action don't match" << endl;
				}
				++amountbreaks;
				Ins::Vertex* vert = solution[d][i];
				if (((ins->t[d].EDT + traveltime[d][i]) - (ins->breakdur) >= ins->breakstart) && ((ins->t[d].EDT + traveltime[d][i]) - ins->breakdur <= ins->breakend))
				{
					breakcheck = true;
					//cout<<"break ok"<<endl;
					//break is scheduled on time
				}
				else
				{
					if (solution[d][i]->index == ins->maxvertices - 1)
					{
						if (ins->t[d].EDT + traveltime[d][i - 1] < ins->breakstart)
						{//route is not long enough to require a break on a regular vertex
							breakcheck = true;
						}
						else
						{
							cout << red << "path" << d << "break wrongly planned at end vertex" << endl;
						}
					}
					else
					{//mistakes to break scheduling found

						if ((ins->t[d].EDT + traveltime[d][i]) - (ins->breakdur) < ins->breakstart)
						{
							cout << red << "path: " << d << "break too early" << endl;
						}
						else
						{
							cout << red << "path: " << d << "break too late" << endl;
						}
					}
				}
			}
		}
		//7. amount of breaks check
		if (amountbreaks == 1)
		{
			//cout << yellow << "break ok" << endl;
		}
		else
		{
			cout << red << "path: " << d << "amount of breaks not ok" << amountbreaks << endl;
		}
		//8. path should start and end at the respective depots
		if ((solution[d][0] == &ins->v[0]) && (solution[d].back() == &ins->v[ins->maxvertices - 1]))
		{
			//cout<<yellow << "start and end ok" << endl;
		}
		else
		{
			cout << red << "start and end vertex not ok" << endl;
		}

	}//end for all d
	for (int i = 1; i < ins->maxvertices; ++i)//enddepot are multiply included
	{
		if (included[i] > 1)
			cout << red << "FAILURE!!! inclusion check fails for vertex:" << i << endl;
	}
	// 6.score check
	if (scorecheck != score)
		cout << red << "score of best solution should be: " << scorecheck << " stored score is: " << score << endl;
}

void Sol::update_traveltime(int tour, int start, int end)//update travel time and max_shift but don't update break
{
	double currenttime = traveltime[tour][start] + ins->t[tour].EDT;
	for (int u = start; u < end-1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = solution[tour][u];
		Ins::Vertex* p = solution[tour][u+1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if (arrivaltime < p->LTW[tour])
		{
			arrivaltime = p->LTW[tour];
		}
		arrivaltime += p->serv + (action[tour][u+1] * ins->breakdur);
		max_shift[tour][u+1] = (traveltime[tour][u+1] + max_shift[tour][u+1]) - (arrivaltime - ins->t[tour].EDT);
		//(old arrival time + old maxshift) - new arrival time = new max_shift
		traveltime[tour][u+1] = arrivaltime - ins->t[tour].EDT;
		currenttime = arrivaltime;
	}//end for
}

void Sol::update_traveltime_break(int tour, int start, int end)//update travel time and maxshift and potentially reschedule break
{
	double reqbreak = true;
	double currenttime = traveltime[tour][start] + ins->t[tour].EDT;
	for (int u = start; u < end - 1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = solution[tour][u];
		Ins::Vertex* p = solution[tour][u + 1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if ((reqbreak) && (arrivaltime >= ins->breakstart))
		{
			action[tour][u+1] = 1;
			breakindex[tour] = u+1;
			reqbreak = false;
		}
		else
		{//erase previously scheduled break
			action[tour][u+1] = 0;
		}
		if (arrivaltime < p->LTW[tour])
		{
			arrivaltime = p->LTW[tour];
		}
		arrivaltime += p->serv + (action[tour][u + 1] * ins->breakdur);
		max_shift[tour][u + 1] = (traveltime[tour][u + 1] + max_shift[tour][u + 1]) - (arrivaltime - ins->t[tour].EDT);
		//(old arrival time + old maxshift) - new arrival time = new max_shift
		traveltime[tour][u + 1] = arrivaltime - ins->t[tour].EDT;
		currenttime = arrivaltime;
	}//end for
	if (reqbreak)
	{
		action[tour].back() = 1;
		traveltime[tour].back() += ins->breakdur;
		breakindex[tour] = (int) solution[tour].size() - 1;
		reqbreak = false;
	}
}

void Sol::update_maxshift(int d, int start, int end, double arrivaltime)
{
	double departuretime = 0;
	Ins::Vertex* y;
	Ins::Vertex* z;
	for (int i = end; i > start ; --i)
	{
		//define the 2 elements
		y = solution[d][i];
		int breaki = action[d][i];
		z = solution[d][i+1];
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
		max_shift[d][i] = departuretime - (traveltime[d][i] + ins->t[d].EDT);
		//reset variable for the calculation of next point
		arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
	}// end for i
}

void Sol::calc_maxshift(int tour)
{
	int size = (int)solution[tour].size();
	double departuretime = 0;
	max_shift[tour].back() = (ins->t[tour].T_max - traveltime[tour].back());
	double arrivaltime = ins->t[tour].LAT - (action[tour].back() * ins->breakdur);//if you break at the end depot subtract breakduration
	Ins::Vertex* y;
	Ins::Vertex* z;
	for (int i = 0; i < size - 2; ++i)// depots don't count
	{
		//define the 2 elements
		y = solution[tour][size - (i + 2)];
		int breaki = action[tour][size - (i + 2)];
		z = solution[tour][size - (i + 1)];
		departuretime = ins->departure_time(y->con[z->index], arrivaltime);
		//TW check on departuretime 
		if (departuretime > y->UTW[tour] + y->serv)
		{
			departuretime = y->UTW[tour] + y->serv;
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
		max_shift[tour][size - (i + 2)] = departuretime - (traveltime[tour][size - (i + 2)] + ins->t[tour].EDT);
		//reset variable for the calculation of next point
		arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
	}// end for i
}

void Sol::calc_maxshift()
{
	for (int d = 0; d < ins->maxtours; ++d)
	{
		int size = (int) solution[d].size();
		double departuretime = 0;
		max_shift[d].back() = (ins->t[d].T_max - traveltime[d].back());
		double arrivaltime = ins->t[d].LAT - (action[d].back() * ins->breakdur);//if you break at the end depot subtract breakduration
		Ins::Vertex* y;
		Ins::Vertex* z;
		for (int i = 0; i < size - 2; ++i)// depots don't count
		{
			//define the 2 elements
			y = solution[d][size - (i + 2)];
			int breaki = action[d][size - (i + 2)];
			z = solution[d][size - (i + 1)];
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
			max_shift[d][size - (i + 2)] = departuretime - (traveltime[d][size - (i + 2)] + ins->t[d].EDT);
			//reset variable for the calculation of next point
			arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
		}// end for i
	}//end for all tours
}

void Sol::insertvertex(int tour, Ins::Vertex* candidate, int position)
{
	available[candidate->index] = false;
	solution[tour].insert(solution[tour].begin() + position + 1, candidate);//insert point y after x
	traveltime[tour].insert(traveltime[tour].begin() + position + 1, 0);//insert temporary value
	action[tour].insert(action[tour].begin() + position + 1, 0);//insert regular visit action change later when necessary
	scores[tour] += candidate->score;// update score of the new solution
	score += candidate->score;// update score of the new solution
	volume[tour] += candidate->volume;
	weight[tour] += candidate->weight;
	if (position < breakindex[tour])
	{
		breakindex[tour] += 1;//due to insertion of 1 vertex the index needs to be incremented with 1
	}
	max_shift[tour].insert(max_shift[tour].begin() + position + 1, 0);
	update_traveltime(tour, position, int(solution[tour].size()));//update travel time and maxshift for all positions after insertion
	double arrivaltime = (traveltime[tour][position + 2] + ins->t[tour].EDT + max_shift[tour][position + 2]) - (solution[tour][position + 2]->serv + action[tour][position + 2] * ins->breakdur);//service time eraftrekken
	update_maxshift(tour, 0, position + 1, arrivaltime);//update maxshift for all positions before insertions
}

void Sol::replacevertex(int tour, Ins::Vertex* candidate, int position)
{
	bool reqbreak = false;
	if (position <= breakindex[tour])
	{
		reqbreak = true;
		for (int vv = 0; vv < solution[tour].size(); ++vv)
		{
			action[tour][vv] = 0;
		}
	}
	Ins::Vertex* old = solution[tour][position];
	available[candidate->index] = false;
	available[old->index] = true;
	solution[tour][position] = candidate;//replace point old with candidate
	max_shift[tour][position] = 0;//dummy value
	score += candidate->score - old->score;// update score of the new solution
	scores[tour] += candidate->score - old->score;// update score of the new solution
	if (reqbreak)
	{
		update_traveltime_break(tour, position - 1, int(solution[tour].size()));
		//if you reposition the break, maxshift has to be recalculated
		calc_maxshift(tour);
	}
	else
	{
		update_traveltime(tour, position - 1, int(solution[tour].size()));//update travel time and maxshift for all positions after replacement
		double arrivaltime = (traveltime[tour][position + 1] + ins->t[tour].EDT + max_shift[tour][position + 1]) - (solution[tour][position + 1]->serv + action[tour][position + 1] * ins->breakdur);//service time eraftrekken
		update_maxshift(tour, 0, position, arrivaltime);//update maxshift for all positions before replacement
	}
	
}

void Sol::removevertex(int tour, int position)
{
	bool reqbreak = false;
	if (position <= breakindex[tour])
	{
		reqbreak = true;
		for (int vv = 0; vv < solution[tour].size(); ++vv)
		{
			action[tour][vv] = 0;
		}
	}
	Ins::Vertex* candidate = solution[tour][position];
	available[candidate->index] = false;
	solution[tour].erase(solution[tour].begin() + position);//insert point y after x
	traveltime[tour].erase(traveltime[tour].begin() + position);//insert temporary value
	action[tour].erase(action[tour].begin() + position);//insert regular visit action change later when necessary
	scores[tour] -= candidate->score;// update score of the new solution
	score -= candidate->score;// update score of the new solution
	volume[tour] -= candidate->volume;
	weight[tour] -= candidate->weight;
	if (reqbreak)
	{
		update_traveltime_break(tour, position-1, int(solution[tour].size()));
	}
	else
	{
		update_traveltime(tour, position-1, int(solution[tour].size()));//update travel time for all after deletion
	}
	calc_maxshift(tour);
}

void Sol::swapvertex(int tour, int i, int j)
{
	Ins::Vertex* remember = solution[tour][i];
	solution[tour][i] = solution[tour][j];
	solution[tour][j] = remember;
	if (i <= breakindex[tour])
	{///break comes after i so might need to be replaced
		update_traveltime_break(tour, i-1, int(solution[tour].size()));
	}
	else
	{
		update_traveltime(tour, i-1, int(solution[tour].size()));
	}
	calc_maxshift(tour);
}

void Sol::optvertices(int tour, int i, int j,bool breakreschedule)
{
	//reverse sequence
	for (int f = 0; f < 1 + (j - i) / 2; ++f)
	{
		Ins::Vertex* temp = solution[tour][j - f];
		solution[tour][j - f] = solution[tour][i + f];
		solution[tour][i + f] = temp;
	}
	if (breakreschedule)
	{
		update_traveltime_break(tour, i - 1, int(solution[tour].size()));
	}
	else
	{
		update_traveltime(tour, i - 1, int(solution[tour].size()));
	}
	calc_maxshift(tour);
}

void Sol::reset() 
{
	for (int t = 0; t < ins->maxtours; ++t)
	{
		solution[t].clear();
		solution[t].push_back(&ins->v[0]);//insert start depot
		traveltime[t].clear();
		traveltime[t].push_back(0.0);//insert first travel time
		action[t].clear();
		action[t].push_back(0);
		max_shift[t].clear();
		max_shift[t].push_back(0);
		scores[t] = 0;
		weight[t] = 0;
		volume[t] = 0;
		breakindex[t] = -1;
	}
	score = 0;
	available.set();//sets all bits to true
	available[ins->v[0].index] = false;
}
