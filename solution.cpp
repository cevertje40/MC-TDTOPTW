#include "solution.h"

Sol::Sol(Ins& ins) :ins(&ins)
{
	tours.resize(ins.maxtours);
	for (int t = 0; t < ins.maxtours; ++t)
	{
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
	shuffle(tourindex.begin(), tourindex.end(), ins.engine);
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
				//cout << "waiting time for: "<<"i"<<i+1<<" , " <<sol.solution[d][i+1]->index<<" <=> "<<waitingtime << endl;
				arrivaltime = current->LTW[d] - (breakcurrent * ins->breakdur);
			}
			arrivaltime += current->serv+breakcurrent*ins->breakdur;
			//cout<<i+1<<" calc traveltime: " << arrivaltime-t[d].EDT << "stored: " << sol.traveltime[d][i+1] << endl;
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
			if (tours[d].deptime[i] + ins->t[d].EDT - tours[d].seq[i]->serv < current->LTW[d])//service time zit al in traveltime
				cout << red << "path: " << d << "FAILURE!!! LTW fail for solutionnr: " << i << " /vertex index: " << current->index << endl;
			if (tours[d].deptime[i] + ins->t[d].EDT - (breakcurrent*ins->breakdur+tours[d].seq[i]->serv) > current->UTW[d])
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
							cout << red << "path" << d << "break wrongly planned at end vertex" << endl;
						}
					}
					else
					{//mistakes to break scheduling found

						if ((ins->t[d].EDT + tours[d].deptime[i]) - (ins->breakdur) < ins->breakstart)
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
	double currenttime = tours[tour].deptime[start] + ins->t[tour].EDT;
	for (int u = start; u < end-1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = tours[tour].seq[u];
		Ins::Vertex* p = tours[tour].seq[u+1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if (arrivaltime + (tours[tour].action[u + 1] * ins->breakdur) < p->LTW[tour])
		{
			arrivaltime = p->LTW[tour]- (tours[tour].action[u + 1] * ins->breakdur);
		}
		arrivaltime += p->serv + (tours[tour].action[u+1] * ins->breakdur);
		tours[tour].max_shift[u+1] = (tours[tour].deptime[u+1] + tours[tour].max_shift[u+1]) - (arrivaltime - ins->t[tour].EDT);
		//(old arrival time + old maxshift) - new arrival time = new max_shift
		tours[tour].deptime[u+1] = arrivaltime - ins->t[tour].EDT;
		currenttime = arrivaltime;
	}//end for
}

void Sol::update_traveltime_break(int tour, int start, int end)//update travel time and maxshift and potentially reschedule break
{
	double reqbreak = true;
	double currenttime = tours[tour].deptime[start] + ins->t[tour].EDT;
	for (int u = start; u < end - 1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = tours[tour].seq[u];
		Ins::Vertex* p = tours[tour].seq[u + 1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if ((reqbreak) && (arrivaltime >= ins->breakstart))
		{
			tours[tour].action[u+1] = 1;
			tours[tour].breakindex = u+1;
			reqbreak = false;
		}
		else
		{//erase previously scheduled break
			tours[tour].action[u+1] = 0;
		}
		if (arrivaltime + (tours[tour].action[u + 1] * ins->breakdur) < p->LTW[tour])
		{
			arrivaltime = p->LTW[tour]- (tours[tour].action[u + 1] * ins->breakdur);
		}
		arrivaltime += p->serv + (tours[tour].action[u + 1] * ins->breakdur);
		tours[tour].max_shift[u + 1] = (tours[tour].deptime[u + 1] + tours[tour].max_shift[u + 1]) - (arrivaltime - ins->t[tour].EDT);
		//(old arrival time + old maxshift) - new arrival time = new max_shift
		tours[tour].deptime[u + 1] = arrivaltime - ins->t[tour].EDT;
		currenttime = arrivaltime;
	}//end for
	if (reqbreak)
	{
		tours[tour].action.back() = 1;
		tours[tour].deptime.back() += ins->breakdur;
		tours[tour].breakindex = (int)tours[tour].seq.size() - 1;
		reqbreak = false;
	}
}

void Sol::update_maxshift(int tour, int start, int end, double arrivaltime)
{
	double departuretime = 0;
	Ins::Vertex* y;
	Ins::Vertex* z;
	for (int i = end; i > start ; --i)
	{
		//define the 2 elements
		y = tours[tour].seq[i];
		int breaki = tours[tour].action[i];
		z = tours[tour].seq[i+1];
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
		tours[tour].max_shift[i] = departuretime - (tours[tour].deptime[i] + ins->t[tour].EDT);
		//reset variable for the calculation of next point
		arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
	}// end for i
}

void Sol::calc_maxshift(Sol::Tour& tour)
{
	int size = (int)tour.seq.size();
	double departuretime = 0;
	tour.max_shift.back() = (ins->t[tour.index].T_max - tour.deptime.back());
	double arrivaltime = ins->t[tour.index].LAT - (tour.action.back() * ins->breakdur);//if you break at the end depot subtract breakduration
	Ins::Vertex* y;
	Ins::Vertex* z;
	for (int i = 0; i < size - 2; ++i)// depots don't count
	{
		//define the 2 elements
		y = tour.seq[size - (i + 2)];
		int breaki = tour.action[size - (i + 2)];
		z = tour.seq[size - (i + 1)];
		departuretime = ins->departure_time(y->con[z->index], arrivaltime);
		//TW check on departuretime 
		if (departuretime > y->UTW[tour.index] + y->serv)
		{
			departuretime = y->UTW[tour.index] + y->serv;
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
		tour.max_shift[size - (i + 2)] = departuretime - (tour.deptime[size - (i + 2)] + ins->t[tour.index].EDT);
		//reset variable for the calculation of next point
		arrivaltime = departuretime - (y->serv + breaki * ins->breakdur);
	}// end for i
}

void Sol::insertvertex(Sol::Tour &tour, Ins::Vertex* candidate, int position)
{
	available[candidate->index] = false;
	tour.seq.insert(tour.seq.begin() + position + 1, candidate);//insert point y after x
	tour.deptime.insert(tour.deptime.begin() + position + 1, 0);//insert temporary value
	tour.action.insert(tour.action.begin() + position + 1, 0);//insert regular visit action change later when necessary
	tour.score += candidate->score;// update score of the new solution
	score += candidate->score;// update score of the new solution
	tour.volume += candidate->volume;
	tour.weight += candidate->weight;
	if (position < tour.breakindex)
	{
		tour.breakindex += 1;//due to insertion of 1 vertex the index needs to be incremented with 1
	}
	tour.max_shift.insert(tour.max_shift.begin() + position + 1, 0);
	update_traveltime(tour.index, position, int(tour.seq.size()));//update travel time and maxshift for all positions after insertion
	double arrivaltime = (tour.deptime[position + 2] + ins->t[tour.index].EDT + tour.max_shift[position + 2]) - (tour.seq[position + 2]->serv + tour.action[position + 2] * ins->breakdur);//service time eraftrekken
	update_maxshift(tour.index, 0, position + 1, arrivaltime);//update maxshift for all positions before insertions
}

void Sol::replacevertex(Sol::Tour &tour, Ins::Vertex* candidate, int position, bool updatebreak)
{
	bool reqbreak = false;
	if (updatebreak)
	{
		if (position <= tour.breakindex)
		{
			reqbreak = true;
			for (int vv = 0; vv < tour.seq.size(); ++vv)
			{
				tour.action[vv] = 0;
			}
		}
	}
	Ins::Vertex* old = tour.seq[position];
	available[candidate->index] = false;
	available[old->index] = true;
	tour.seq[position] = candidate;//replace point old with candidate
	tour.max_shift[position] = 0;//dummy value
	tour.volume -= old->volume - candidate->volume;
	tour.weight -= old->weight - candidate->weight;
	score += candidate->score - old->score;// update score of the new solution
	tour.score += candidate->score - old->score;// update score of the new solution
	if (reqbreak)
	{
		update_traveltime_break(tour.index, position - 1, int(tour.seq.size()));
		//if you reposition the break, maxshift has to be recalculated
		calc_maxshift(tour);
	}
	else
	{
		update_traveltime(tour.index, position - 1, int(tour.seq.size()));//update travel time and maxshift for all positions after replacement
		double arrivaltime = (tour.deptime[position + 1] + ins->t[tour.index].EDT + tour.max_shift[position + 1]) - (tour.seq[position + 1]->serv + tour.action[position + 1] * ins->breakdur);//service time eraftrekken
		update_maxshift(tour.index, 0, position, arrivaltime);//update maxshift for all positions before replacement
	}
	
}

void Sol::removevertex(Sol::Tour& tour, int position)
{
	bool reqbreak = false;
	if (position <= tour.breakindex)
	{
		reqbreak = true;
		for (int vv = 0; vv < tour.seq.size(); ++vv)
		{
			tour.action[vv] = 0;
		}
	}
	Ins::Vertex* candidate = tour.seq[position];
	available[candidate->index] = false;
	tour.seq.erase(tour.seq.begin() + position);//insert point y after x
	tour.deptime.erase(tour.deptime.begin() + position);//insert temporary value
	tour.action.erase(tour.action.begin() + position);//insert regular visit action change later when necessary
	tour.score -= candidate->score;// update score of the new solution
	score -= candidate->score;// update score of the new solution
	tour.volume -= candidate->volume;
	tour.weight -= candidate->weight;
	if (reqbreak)
	{
		update_traveltime_break(tour.index, position-1, int(tour.seq.size()));
	}
	else
	{
		update_traveltime(tour.index, position-1, int(tour.seq.size()));//update travel time for all after deletion
	}
	calc_maxshift(tour);
}

void Sol::swapvertex(Sol::Tour &tour, int i, int j)
{
	Ins::Vertex* remember = tour.seq[i];
	tour.seq[i] = tour.seq[j];
	tour.seq[j] = remember;
	bool reqbreak = false;
	if (i <= tour.breakindex)
	{
		reqbreak = true;
		for (int vv = i; vv < tour.seq.size(); ++vv)
		{
			tour.action[vv] = 0;
		}
	}
	if (reqbreak)
	{///break comes after i so might need to be replaced
		update_traveltime_break(tour.index, i-1, int(tour.seq.size()));
	}
	else
	{
		update_traveltime(tour.index, i-1, int(tour.seq.size()));
	}
	calc_maxshift(tour);
}

void Sol::optvertices(Sol::Tour &tour, int i, int j)
{
	//reverse sequence
	for (int f = 0; f < 1 + (j - i) / 2; ++f)
	{
		Ins::Vertex* temp = tour.seq[j - f];
		tour.seq[j - f] = tour.seq[i + f];
		tour.seq[i + f] = temp;
	}
	bool reqbreak = false;
	if (i <= tour.breakindex)
	{
		reqbreak = true;
		for (int vv = i; vv < tour.seq.size(); ++vv)
		{
			tour.action[vv] = 0;
		}
	}
	if (reqbreak)
	{
		update_traveltime_break(tour.index, i - 1, int(tour.seq.size()));
	}
	else
	{
		update_traveltime(tour.index, i - 1, int(tour.seq.size()));
	}
	calc_maxshift(tour);
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
