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
	}
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
