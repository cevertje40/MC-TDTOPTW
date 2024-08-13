#include "solution.h"

Sol::Sol(Ins& ins):ins(&ins)
{
	tours.resize(ins.maxtours);
	for (int t = 0; t < ins.maxtours; ++t)
	{
		tours[t].ins = &ins;
		tours[t].seq.reserve(ins.maxvertices);//reserves memory
		tours[t].seq.push_back(&ins.v[0]);//insert start depot
		tours[t].deptime.reserve(ins.maxvertices);//reserve memory
		tours[t].deptime.push_back(0);//insert first travel time
		tours[t].action.reserve(ins.maxvertices);
		tours[t].action.push_back(0);
		tours[t].max_shift.reserve(ins.maxvertices);
		tours[t].max_shift.push_back(0);
		tours[t].score = 0;
		tours[t].weight = 0.0;
		tours[t].volume = 0.0;
		tours[t].breakindex = -1;
		tours[t].index = t;
		tourindex.push_back(t);
	}
	shuffle(tourindex.begin(), tourindex.end(),engine);
	score = 0;
	available = boost::dynamic_bitset<>(ins.maxvertices);
	available.set();//sets all bits to true
	available[ins.v[0].index] = false;
}

ostream& operator<<(ostream& output, Sol& sol)
{
	const char sep = ' ';
	for (int d = 0; d < (int)sol.tours.size(); ++d)
	{
		output << "vehicle: ";
		output << left << setw(2) << setfill(sep) << d;
		output << " score: ";
		output << left << setw(3) << setfill(sep) << sol.tours[d].score;
		output << " weight: ";
		output << left << setw(5) << setfill(sep) << sol.tours[d].weight << "/" << sol.ins->t[d].W_max;
		output << " volume: ";
		output << left << setw(5) << setfill(sep) << sol.tours[d].volume << "/" << sol.ins->t[d].V_max;
		output << " traveltime: ";
		output << left << setw(5) << setfill(sep) << sol.tours[d].deptime.back() << "/" << sol.ins->t[d].T_max << "\n";
	}
	for (int d = 0; d < (int)sol.tours.size(); ++d)
	{
		int end = (int)sol.tours[d].seq.size();
		output << left << setw(4) << setfill(sep) << "d" << d << "\n";
		output << left << setw(4) << setfill(sep) << "i";
		output << left << setw(4) << setfill(sep) << "vi";
		output << left << setw(7) << setfill(sep) << "sco";
		output << left << setw(7) << setfill(sep) << "vol";
		output << left << setw(7) << setfill(sep) << "wei";
		output << left << setw(8) << setfill(sep) << "dep";
		output << left << setw(8) << setfill(sep) << "shift" << "\n";
		for (int i = 0; i < end; ++i)
		{
			output << left << setw(4) << setfill(sep) << i;
			output << left << setw(4) << setfill(sep) << sol.tours[d].seq[i]->index;
			output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->score;
			output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->volume;
			output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->weight;
			output << left << setw(8) << setfill(sep) << sol.tours[d].deptime[i] << " [" << sol.tours[d].seq[i]->LTW[d] - sol.ins->t[d].EDT << " ; " << sol.tours[d].seq[i]->UTW[d] - sol.ins->t[d].EDT << "]";
			output << left << setw(8) << setfill(sep) << sol.tours[d].max_shift[i] << "\n";
		}
	}
	output << "total score : " << sol.score;
	return output;
}//end output operator

void Sol::check()
{
	int scorecheck = 0;
	vector<int> included(ins->maxvertices, 0);
	for (int d=0;d<ins->maxtours;++d)
	{
		//cout << "tour: " << d << endl << endl;
		double weightcheck = 0;
		double volumecheck = 0;
		//1. travel time check
		double currenttime = ins->t[d].EDT + tours[d].deptime[0];
		double length = 0;
		int end = (int)tours[d].seq.size() - 1;
		for (int i = 0; i < end; ++i)
		{
			++included[tours[d].seq[i]->index];
			scorecheck += tours[d].seq[i]->score;
			weightcheck += tours[d].seq[i]->weight;
			volumecheck += tours[d].seq[i]->volume;
			Ins::Vertex* last = tours[d].seq[i];
			Ins::Vertex* current = tours[d].seq[i + 1];
			int breakcurrent = tours[d].action[i + 1];
			double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
			double waitingtime = 0;
			if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[d])
			{
				waitingtime = current->LTW[d] - (arrivaltime + breakcurrent * (ins->breakdur));
				//cout << "waiting time for: "<<"i"<<i+1<<" , " <<tours[d].seq[i + 1]->index << " <=> " << waitingtime << endl;
				arrivaltime = current->LTW[d] - (breakcurrent * ins->breakdur);
			}
			arrivaltime += current->serv+breakcurrent*ins->breakdur;
			//cout<<i+1<<" calc traveltime: " << arrivaltime-ins->t[d].EDT << " stored: " << tours[d].deptime[i + 1] << endl;
			currenttime = arrivaltime;
		}//end for i
		length = currenttime - ins->t[d].EDT;
		if (abs(length-tours[d].deptime.back())>0.01)
			cout << red << "path: " << d << "new calculated length: " << length << " stored length: " << tours[d].deptime.back() << "max length" << ins->t[d].T_max << endl;
		//3. TW check
		for (int i = 0; i <= end; ++i)//utw van end depot ook checken
		{
			Ins::Vertex* current = tours[d].seq[i];
			int breakcurrent = tours[d].action[i];
			if (tours[d].deptime[i] + ins->t[d].EDT - tours[d].seq[i]->serv+0.01 < current->LTW[d])//service time zit al in traveltime
				cout << red << "path: " << d << "FAILURE!!! LTW fail for solutionnr: " << i << " /vertex index: " << current->index << endl;
			if (((tours[d].deptime[i] + ins->t[d].EDT) - tours[d].seq[i]->serv) -0.01> current->UTW[d])
			{
				cout << red << "path: " << d << "FAILURE!!! UTW fail for solutionnr: " << i << " /vertex index: " << current->index << endl;
			}
		}//end for i
		//4. weight check
		if (abs(weightcheck-tours[d].weight)>0.01)
			cout << red << "path: " << d << "new calculated weight" << weightcheck << "stored weight: " << tours[d].weight << "max: " << ins->t[d].W_max << endl;
		//5. volume check
		if (abs(volumecheck-tours[d].volume)>0.01)
			cout << yellow << "path: " << d << "new calculated volume" << volumecheck << "stored volume " << tours[d].volume << "max: " << ins->t[d].V_max << endl;
		//6. break timing check
		bool breakcheck = false;
		int amountbreaks = 0;
		for (int i = 0; i <= end; ++i)//break kan op enddepot zitten
		{
			if (tours[d].action[i] == 1)
			{
				if (tours[d].breakindex != i)
				{
					cout << red << "path: " << d << "breakindex and sol action don't match" << endl;
				}
				++amountbreaks;
				Ins::Vertex* vert = tours[d].seq[i];
				//check if break starts within allowed time frame, break is taken before service as break can be taken before opening of vertex
				if (((ins->t[d].EDT + tours[d].deptime[i]-tours[d].seq[i]->serv) - (ins->breakdur) >= ins->breakstart) && ((ins->t[d].EDT + tours[d].deptime[i] - tours[d].seq[i]->serv) - ins->breakdur <= ins->breakend))
				{
					breakcheck = true;
					//cout<<"break ok"<<endl;
					//break is scheduled on time
				}
				else
				{
					if (tours[d].seq[i]->index == ins->maxvertices - 1)
					{
						if (ins->t[d].EDT + tours[d].deptime[i - 1] < ins->breakstart)
						{//route is not long enough to require a break on a regular vertex
							breakcheck = true;
						}
						else
						{
							cout << red << "path" << d << " break wrongly planned at end vertex" << endl;
						}
					}
					else
					{//mistakes to break scheduling found

						if ((ins->t[d].EDT + tours[d].deptime[i] - tours[d].seq[i]->serv) - (ins->breakdur) < ins->breakstart)
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
		if ((tours[d].seq[0] == &ins->v[0]) && (tours[d].seq.back() == &ins->v[ins->maxvertices - 1]))
		{
			//cout<<yellow << "start and end ok" << endl;
		}
		else
		{
			cout << red << "start and end vertex not ok" << endl;
		}
		// 9.max_shift check
		vector<double> max_shiftcheck(tours[d].max_shift.size(), 0);
		int size = (int)tours[d].seq.size();
		double departuretime = 0;
		max_shiftcheck.back() = (ins->t[d].T_max - tours[d].deptime.back());
		double arrivaltime = ins->t[d].LAT - (tours[d].action.back() * ins->breakdur);//if you break at the end depot subtract breakduration
		if (tours[d].action.back() == 1)//break op enddepot
		{
			if (ins->t[d].EDT + ins->t[d].T_max > ins->breakend + ins->breakdur)
			{
				//cout<<"path: "<<d<< " bij calc maxshift break op enddepot verhindert een maxshift: " << endl;
				max_shiftcheck.back() = (ins->breakend + ins->breakdur) - (tours[d].deptime.back() + ins->t[d].EDT);
				arrivaltime = ins->breakend;//zoals hieronder service of enkel break in dit geval ervan aftrekken
			}
		}
		Ins::Vertex* y;
		Ins::Vertex* z;
		for (int i = 0; i < size - 2; ++i)// depots don't count
		{
			//define the 2 vertices
			y = tours[d].seq[size - (i + 2)];
			int breaki = tours[d].action[size - (i + 2)];
			z = tours[d].seq[size - (i + 1)];
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
					departuretime = ins->breakend + ins->breakdur;
				}
			}
			//store result
			max_shiftcheck[size - (i + 2)] = departuretime - (tours[d].deptime[size - (i + 2)] + ins->t[d].EDT);
			//reset variable for the calculation of next point
			arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
		}// end for i
		for (int i = 0; i <= end; ++i)//break kan op enddepot zitten
		{
			if (tours[d].max_shift[i] != max_shiftcheck[i])
			{
				cout << red << "error in max_shift for position: " << i << endl;
			}
		}
	}//end for all tours

	//10. inclusion & availability check
	for (int i = 1; i < ins->maxvertices-1; ++i)//depot can be included multiple times
	{
		if (included[i] > 1)//regular vertex included more than once
		{
			cout << red << "FAILURE!!! inclusion check fails for vertex:" << i << endl;
		}
		if (included[i] == 1)
		{
			if (available[i] == true)
			{
				cout << red << "error in availability bitset for vertex: " << i << endl;
			}
		}
		if (included[i] == 0)
		{
			if (available[i] == false)
			{
				cout << red << "error in availability bitset for vertex: " << i << endl;
			}
		}
	}
	// 11. total score check
	if (scorecheck != score)
	{
		cout << red << "score of best solution should be: " << scorecheck << " stored score is: " << score << endl;
	}
}

void Tour::update_traveltime(int start, int end)//update travel time and max_shift but don't update break
{
	double currenttime = deptime[start] + ins->t[index].EDT;
	for (int u = start; u < end-1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u+1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if (arrivaltime + (action[u + 1] * ins->breakdur) < p->LTW[index])
		{
			arrivaltime = p->LTW[index]- (action[u + 1] * ins->breakdur);
		}
		arrivaltime += p->serv + (action[u+1] * ins->breakdur);
		max_shift[u+1] = (deptime[u+1] + max_shift[u+1]) - (arrivaltime - ins->t[index].EDT);
		//(old arrival time + old maxshift) - new arrival time = new max_shift
		deptime[u+1] = arrivaltime - ins->t[index].EDT;
		currenttime = arrivaltime;
	}//end for
}

void Tour::update_traveltime_break(int start, int end)//update travel time and maxshift and potentially reschedule break
{
	//evaluate tour without break
	for (int vv = 0; vv < action.size(); ++vv)
	{
		action[vv] = 0;
	}
	double currenttime = ins->t[index].EDT;
	int earliestbreakindex = (int) seq.size()-1;//default to end depot
	bool earliestbreaktofind = true;
	vector<double>waitingtime(seq.size(),0.0);
	for (int u = 0; u < end - 1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if ((earliestbreaktofind) && (max(p->LTW[index] - ins->breakdur, arrivaltime) >= ins->breakstart))//arrival time need to be higher than breakstart not dep-serv
		{//determine earliest feasible point to insert a break
			earliestbreakindex = u+1;
			earliestbreaktofind = false;
		}
		if (arrivaltime < p->LTW[index])
		{
			waitingtime[u+1] = p->LTW[index] - arrivaltime;
			arrivaltime = p->LTW[index];

		}
		arrivaltime += p->serv;
		deptime[u + 1] = arrivaltime - ins->t[index].EDT;
		currenttime = arrivaltime;
	}//end for
	//calc max_shift to check where there is room to insert a break
	calc_maxshift();
	//schedule the break 
	bool breakindexnotfound = true;
	for (int u = earliestbreakindex; u <= end-1; ++u)//end depot included
	{
		if (max_shift[u]+waitingtime[u] >= ins->breakdur)
		{
			action[u] = 1;
			breakindex = u;
			deptime[u] += ins->breakdur;
			breakindexnotfound = false;
			break;
		}
	}
	if (breakindexnotfound)
	{//put break on enddepot if the route is too short
		action.back() = 1;
		
	}

	if (breakindex == seq.size() - 1)
	{
		double arrivaltime;
		if (ins->t[index].EDT + ins->t[index].T_max > ins->breakend + ins->breakdur)
		{
			max_shift.back() = (ins->breakend + ins->breakdur) - (deptime.back() + ins->t[index].EDT);
			arrivaltime = ins->breakend;
		}
		update_maxshift(0, breakindex - 1, arrivaltime);
	}
	else
	{
		//update traveltime and max_shift for vertices after breakindex-1
		currenttime = deptime[breakindex - 1] + ins->t[index].EDT;
		for (int u = breakindex - 1; u < end - 1; ++u)
		{
			Ins::Vertex* o = seq[u];
			Ins::Vertex* p = seq[u + 1];
			//calc travel time from van o to p
			double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
			//LTW check
			if (arrivaltime + (action[u + 1] * ins->breakdur) < p->LTW[index])
			{
				arrivaltime = p->LTW[index] - (action[u + 1] * ins->breakdur);
			}
			arrivaltime += p->serv + (action[u + 1] * ins->breakdur);
			max_shift[u + 1] = (deptime[u + 1] + max_shift[u + 1]) - (arrivaltime - ins->t[index].EDT);
			deptime[u + 1] = arrivaltime - ins->t[index].EDT;
			currenttime = arrivaltime;
		}
		//update max_shift for vertices before breakindex
		double arrivaltime = (deptime[breakindex + 1] + ins->t[index].EDT + max_shift[breakindex + 1]) - (seq[breakindex + 1]->serv + action[breakindex + 1] * ins->breakdur);//service time eraftrekken
		update_maxshift(0, breakindex, arrivaltime);//update maxshift for all positions before breakindex
	}
}

void Tour::update_maxshift(int start, int end, double arrivaltime)
{
	double departuretime = 0;
	Ins::Vertex* y;
	Ins::Vertex* z;
	for (int i = end; i > start ; --i)
	{
		//define the 2 elements
		y = seq[i];
		int breaki = action[i];
		z = seq[i+1];
		departuretime = ins->departure_time(y->con[z->index], arrivaltime);
		//TW check on departuretime 
		if (departuretime > y->UTW[index] + y->serv)
		{
			departuretime = y->UTW[index] + y->serv;
		}
		//break may limit the maximum allowable shift
		if (breaki == 1)
		{
			if (departuretime > ins->breakend + ins->breakdur)
			{
				//cout << "bij calc maxshift break verhindert een maxshift: " << departuretime << "<=>" << ins->breakend + ins->breakdur << endl;
				departuretime = ins->breakend + ins->breakdur;
			}
		}
		//store result
		max_shift[i] = departuretime - (deptime[i] + ins->t[index].EDT);
		//reset variable for the calculation of next point
		arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
	}// end for i
}//end update max_shift

void Tour::calc_maxshift()
{
	int size = (int) seq.size();
	double departuretime = 0;
	max_shift.back() = (ins->t[index].T_max - deptime.back());
	double arrivaltime = ins->t[index].LAT - (action.back() * ins->breakdur);//if you break at the end depot subtract breakduration
	if (action.back() == 1)//break op enddepot
	{
		if (ins->t[index].EDT+ins->t[index].T_max > ins->breakend + ins->breakdur)
		{
			//cout<<"path: "<<d<< " bij calc maxshift break op enddepot verhindert een maxshift: " << endl;
			max_shift.back() = (ins->breakend + ins->breakdur) - (deptime.back() + ins->t[index].EDT);
			arrivaltime = ins->breakend;//zoals hieronder service of enkel break in dit geval ervan aftrekken
		}
	}
	Ins::Vertex* y;
	Ins::Vertex* z;
	for (int i = 0; i < size - 2; ++i)// depots don't count
	{
		//define the 2 elements
		y = seq[size - (i + 2)];
		int breaki = action[size - (i + 2)];
		z = seq[size - (i + 1)];
		departuretime = ins->departure_time(y->con[z->index], arrivaltime);
		//TW check on departuretime 
		if (departuretime > y->UTW[index] + y->serv)
		{
			departuretime = y->UTW[index] + y->serv;
		}
		//break may limit the maximum allowable shift
		if (breaki == 1)
		{
			if (departuretime > ins->breakend + ins->breakdur)
			{
				//cout << "bij calc maxshift break verhindert een maxshift: " << departuretime << "<=>" << ins->breakend + ins->breakdur << endl;
				departuretime = ins->breakend + ins->breakdur;
			}
		}
		//store result
		max_shift[size - (i + 2)] = departuretime - (deptime[size - (i + 2)] + ins->t[index].EDT);
		//reset variable for the calculation of next point
		arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
	}// end for i
}//end calc_max_shift

void Tour::insert_vertex(Ins::Vertex* candidate, int position)
{
	seq.insert(seq.begin() + position + 1, candidate);//insert point y after x
	deptime.insert(deptime.begin() + position + 1, 0);//insert temporary value
	action.insert(action.begin() + position + 1, 0);//insert regular visit action change later when necessary
	score += candidate->score;// update score of the new solution
	volume += candidate->volume;
	weight += candidate->weight;
	if (position < breakindex)
	{
		breakindex += 1;//due to insertion of 1 vertex the index needs to be incremented with 1
	}
	max_shift.insert(max_shift.begin() + position + 1, 0);
	update_traveltime(position, int(seq.size()));//update travel time and maxshift for all positions after insertion
	double arrivaltime = (deptime[position + 2] + ins->t[index].EDT + max_shift[position + 2]) - (seq[position + 2]->serv + action[position + 2] * ins->breakdur);//service time eraftrekken
	update_maxshift(0, position + 1, arrivaltime);//update maxshift for all positions before insertions
}

void Tour::remove_vertex(int position)
{
	bool reqbreak = false;
	if (position <= breakindex)
	{
		reqbreak = true;
		for (int vv = 0; vv < seq.size(); ++vv)
		{
			action[vv] = 0;
		}
	}
	Ins::Vertex* candidate = seq[position];
	seq.erase(seq.begin() + position);//insert point y after x
	deptime.erase(deptime.begin() + position);//insert temporary value
	max_shift.erase(max_shift.begin() + position);
	action.erase(action.begin() + position);//insert regular visit action change later when necessary
	score -= candidate->score;// update score of the tour
	volume -= candidate->volume;
	weight -= candidate->weight;
	if (reqbreak)
	{
		update_traveltime_break(position - 1, int(seq.size()));
	}
	else
	{
		update_traveltime(position - 1, int(seq.size()));//update travel time for all after deletion
	}
	calc_maxshift();
}//end remove_vertex tour version

void Tour::remove_vertices(int position1, int position2)
{
	Ins::Vertex* candidate1 = seq[position1];
	Ins::Vertex* candidate2 = seq[position2];
	score -= candidate1->score;// update score of the tour
	score -= candidate2->score;// update score of the tour
	volume -= candidate1->volume;
	weight -= candidate1->weight;
	volume -= candidate2->volume;
	weight -= candidate2->weight;
	seq.erase(seq.begin() + position1);//insert point y after x
	deptime.erase(deptime.begin() + position1);//insert temporary value
	max_shift.erase(max_shift.begin() + position1);
	action.erase(action.begin() + position1);//insert regular visit action change later when necessary
	if (position1 < position2)
	{
		position2 -= 1;
	}
	seq.erase(seq.begin() + position2);//insert point y after x
	deptime.erase(deptime.begin() + position2);//insert temporary value
	max_shift.erase(max_shift.begin() + position2);
	action.erase(action.begin() + position2);//insert regular visit action change later when necessary
	update_traveltime_break(position1 - 1, int(seq.size()));
}

void Tour::replace_vertex(Ins::Vertex* candidate, int position)
{
	bool reqbreak = false;
	if (position <= breakindex)
	{
		reqbreak = true;
	}
	Ins::Vertex* old = seq[position];
	seq[position] = candidate;//replace point old with candidate
	max_shift[position] = 0;//dummy value
	volume -= old->volume - candidate->volume;
	weight -= old->weight - candidate->weight;
	score += candidate->score - old->score;// update score of the new solution
	if (reqbreak)
	{
		update_traveltime_break(position - 1, int(seq.size()));
		//if you reposition the break, maxshift has to be recalculated
		calc_maxshift();
	}
	else
	{
		update_traveltime(position - 1, int(seq.size()));//update travel time and maxshift for all positions after replacement
		double arrivaltime = (deptime[position + 1] + ins->t[index].EDT + max_shift[position + 1]) - (seq[position + 1]->serv + action[position + 1] * ins->breakdur);//service time eraftrekken
		update_maxshift(0, position, arrivaltime);//update maxshift for all positions before replacement
	}
}

void Tour::opt_vertices(int i, int j)
{
	//reverse sequence
	for (int f = 0; f < 1 + (j - i) / 2; ++f)
	{
		Ins::Vertex* temp = seq[j - f];
		seq[j - f] = seq[i + f];
		seq[i + f] = temp;
	}
	//determine wether break is required
	bool reqbreak = false;
	if (i <= breakindex)
	{
		reqbreak = true;
	}
	if (reqbreak)
	{
		update_traveltime_break(i - 1, int(seq.size()));
	}
	else
	{
		update_traveltime(i - 1, int(seq.size()));
	}
	calc_maxshift();
}

void Tour::swap_vertices(int i, int j)
{
	Ins::Vertex* remember = seq[i];
	seq[i] = seq[j];
	seq[j] = remember;
	bool reqbreak = false;
	if (i <= breakindex)
	{
		reqbreak = true;
	}
	if (reqbreak)
	{///break comes after i so might need to be replaced
		update_traveltime_break(i - 1, int(seq.size()));
	}
	else
	{
		update_traveltime(i - 1, int(seq.size()));
	}
	calc_maxshift();
}


void Tour::check()
{
	int scorecheck = 0;
	//cout << "tour: " << d << endl << endl;
	double weightcheck = 0;
	double volumecheck = 0;
	//1. travel time check
	double currenttime = ins->t[index].EDT + deptime[0];
	double length = 0;
	int end = (int)seq.size() - 1;
	for (int i = 0; i < end; ++i)
	{
		scorecheck += seq[i]->score;
		weightcheck += seq[i]->weight;
		volumecheck += seq[i]->volume;
		Ins::Vertex* last = seq[i];
		Ins::Vertex* current = seq[i + 1];
		int breakcurrent = action[i + 1];
		double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);
		double waitingtime = 0;
		if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[index])
		{
			waitingtime = current->LTW[index] - (arrivaltime + breakcurrent * (ins->breakdur));
			//cout << "waiting time for: "<<"i"<<i+1<<" , " <<tours[d].seq[i + 1]->index << " <=> " << waitingtime << endl;
			arrivaltime = current->LTW[index] - (breakcurrent * ins->breakdur);
		}
		arrivaltime += current->serv + breakcurrent * ins->breakdur;
		//cout<<i+1<<" calc traveltime: " << arrivaltime-ins->t[d].EDT << " stored: " << tours[d].deptime[i + 1] << endl;
		currenttime = arrivaltime;
	}//end for i
	length = currenttime - ins->t[index].EDT;
	if (abs(length - deptime.back()) > 0.01)
		cout << red << "tour: " << index << "new calculated length: " << length << " stored length: " << deptime.back() << "max length" << ins->t[index].T_max << endl;
	//3. TW check
	for (int i = 0; i <= end; ++i)//utw van end depot ook checken
	{
		Ins::Vertex* current = seq[i];
		int breakcurrent = action[i];
		if (deptime[i] + ins->t[index].EDT - seq[i]->serv + 0.01 < current->LTW[index])//service time zit al in traveltime
			cout << red << "tour: " << index << "FAILURE!!! LTW fail for solutionnr: " << i << " /vertex index: " << current->index << endl;
		if (((deptime[i] + ins->t[index].EDT) - seq[i]->serv) - 0.01 > current->UTW[index])
		{
			cout << red << "tour: " << index << "FAILURE!!! UTW fail for solutionnr: " << i << " /vertex index: " << current->index << endl;
		}
	}//end for i
	//4. weight check
	if (abs(weightcheck - weight) > 0.01)
		cout << red << "tour: " << index << "new calculated weight" << weightcheck << "stored weight: " << weight << "max: " << ins->t[index].W_max << endl;
	//5. volume check
	if (abs(volumecheck - volume) > 0.01)
		cout << yellow << "tour: " << index << "new calculated volume" << volumecheck << "stored volume " << volume << "max: " << ins->t[index].V_max << endl;
	//6. break timing check
	bool breakcheck = false;
	int amountbreaks = 0;
	for (int i = 0; i <= end; ++i)//break kan op enddepot zitten
	{
		if (action[i] == 1)
		{
			if (breakindex != i)
			{
				cout << red << "path: " << index << "breakindex and sol action don't match" << endl;
			}
			++amountbreaks;
			Ins::Vertex* vert = seq[i];
			//check if break starts within allowed time frame, break is taken before service as break can be taken before opening of vertex
			if (((ins->t[index].EDT + deptime[i] - seq[i]->serv) - (ins->breakdur) >= ins->breakstart) && ((ins->t[index].EDT + deptime[i] - seq[i]->serv) - ins->breakdur <= ins->breakend))
			{
				breakcheck = true;
				//cout<<"break ok"<<endl;
				//break is scheduled on time
			}
			else
			{
				if (seq[i]->index == ins->maxvertices - 1)
				{
					if (ins->t[index].EDT + deptime[i - 1] < ins->breakstart)
					{//route is not long enough to require a break on a regular vertex
						breakcheck = true;
					}
					else
					{
						cout << red << "tour: " << index << " break wrongly planned at end vertex" << endl;
					}
				}
				else
				{//mistakes to break scheduling found

					if ((ins->t[index].EDT + deptime[i] - seq[i]->serv) - (ins->breakdur) < ins->breakstart)
					{
						cout << red << "tour: " << index << "break too early" << endl;
					}
					else
					{
						cout << red << "tour: " << index << "break too late" << endl;
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
		cout << red << "tour: " << index << "amount of breaks not ok" << amountbreaks << endl;
	}
	//8. tour should start and end at the respective depots
	if ((seq[0] == &ins->v[0]) && (seq.back() == &ins->v[ins->maxvertices - 1]))
	{
		//cout<<yellow << "start and end ok" << endl;
	}
	else
	{
		cout << red<<"tour: " << index << "start and end vertex not ok" << endl;
	}
	// 9.max_shift check
	vector<double> max_shiftcheck(max_shift.size(), 0);
	int size = (int)seq.size();
	double departuretime = 0;
	max_shiftcheck.back() = (ins->t[index].T_max - deptime.back());
	double arrivaltime = ins->t[index].LAT - (action.back() * ins->breakdur);//if you break at the end depot subtract breakduration
	if (action.back() == 1)//break op enddepot
	{
		if (ins->t[index].EDT + ins->t[index].T_max > ins->breakend + ins->breakdur)
		{
			//cout<<"path: "<<d<< " bij calc maxshift break op enddepot verhindert een maxshift: " << endl;
			max_shiftcheck.back() = (ins->breakend + ins->breakdur) - (deptime.back() + ins->t[index].EDT);
			arrivaltime = ins->breakend;//zoals hieronder service of enkel break in dit geval ervan aftrekken
		}
	}
	Ins::Vertex* y;
	Ins::Vertex* z;
	for (int i = 0; i < size - 2; ++i)// depots don't count
	{
		//define the 2 elements
		y = seq[size - (i + 2)];
		int breaki = action[size - (i + 2)];
		z = seq[size - (i + 1)];
		departuretime = ins->departure_time(y->con[z->index], arrivaltime);
		//TW check on departuretime 
		if (departuretime > y->UTW[index] + y->serv)
		{
			departuretime = y->UTW[index] + y->serv;
		}
		//break may limit the maximum allowable shift
		if (breaki == 1)
		{
			if (departuretime > ins->breakend + ins->breakdur)
			{
				//cout << "bij calc maxshift break verhindert een maxshift: " << departuretime << "<=>" << ins->breakend + ins->breakdur << endl;
				departuretime = ins->breakend + ins->breakdur;
			}
		}
		//store result
		max_shiftcheck[size - (i + 2)] = departuretime - (deptime[size - (i + 2)] + ins->t[index].EDT);
		//reset variable for the calculation of next point
		arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
	}// end for i
	for (int i = 0; i <= end; ++i)
	{
		if (max_shift[i] != max_shiftcheck[i])
		{
			cout << red <<"tour: "<< index << "error in max_shift for position: " << i <<" new max_shift: "<<max_shiftcheck[i]<<" stored max_shift: "<<max_shift[i] << endl;
		}
	}
}//end tour check


void Sol::insert_vertex(Sol::Tour &tour, Ins::Vertex* candidate, int position)
{
	available[candidate->index] = false;
	score += candidate->score;// update score of the new solution
	tour.insert_vertex(candidate, position);
}

void Sol::replace_vertex(Sol::Tour &tour, Ins::Vertex* candidate, int position)
{
	Ins::Vertex* old = tour.seq[position];
	available[candidate->index] = false;
	available[old->index] = true;
	score += candidate->score - old->score;// update score of the new solution
	tour.replace_vertex(candidate, position);
}

void Sol::remove_vertex(Sol::Tour& tour, int position)
{
	Ins::Vertex* candidate = tour.seq[position];
	available[candidate->index] = true;//make vertex back available
	score -= candidate->score;// update score of the new solution
	tour.remove_vertex(position);
}

void Sol::remove_vertices(Tour& tour, int position1, int position2)
{
	Ins::Vertex* candidate1 = tour.seq[position1];
	Ins::Vertex* candidate2 = tour.seq[position2];
	available[candidate1->index] = true;//make vertex back available
	score -= candidate1->score;// update score of the new solution
	available[candidate2->index] = true;//make vertex back available
	score -= candidate2->score;// update score of the new solution
	tour.remove_vertices(position1,position2);
}

void Sol::reset() 
{
	for (int tour = 0; tour < ins->maxtours; ++tour)
	{
		tours[tour].seq.clear();
		tours[tour].seq.push_back(&ins->v[0]);//insert start depot
		tours[tour].deptime.clear();
		tours[tour].deptime.push_back(0.0);//insert first travel time
		tours[tour].action.clear();
		tours[tour].action.push_back(0);
		tours[tour].max_shift.clear();
		tours[tour].max_shift.push_back(0);
		tours[tour].score = 0;
		tours[tour].weight = 0;
		tours[tour].volume = 0;
		tours[tour].breakindex = -1;
	}
	score = 0;
	available.set();//sets all bits to true
	available[ins->v[0].index] = false;
}
