#pragma once
#include "instance.h"



class Sol
{
public:
	class Tour
	{
		public:
		int index;
		vector <Ins::Vertex*> seq;//sequence of vertex pointers
		vector<double> deptime;//departuretime-EDT at each vertex
		vector<double> max_shift;//local evaluation metric, maximum amount of time each vertex can be shifted forward in time
		vector<int> action;//0 visit, 1 break and visit
		int score;//total score of tour
		double weight;//weight per tour
		double volume;//volume per tour
		int breakindex;//position of break in tour
	};
	vector<Tour>tours;//solution consist of collection of tours
	vector<int>tourindex;//random tour index
	boost::dynamic_bitset<> available;// bitset that states for every vertex if it is still available for inclusion
	int score;//sum of all tour scores
	Ins* ins;//pointer to instance object
	//methods
	Sol(){}
	~Sol(){}
	Sol(Ins& ins);
	void reset();
	void check();
	void update_traveltime(int tour, int start, int end);
	void update_traveltime_break(int tour, int start, int end);
	void update_maxshift(int tour, int start, int end, double arrivaltime);//partial update within a tour
	void calc_maxshift(int tour);//for specific tour of solution
	void calc_maxshift(Sol::Tour &tour);//for specific tour of solution
	void calc_maxshift();//for all tours for whole solution
	void insertvertex(int tour, Ins::Vertex* candidate, int position);
	void replacevertex(int tour, Ins::Vertex* candidate, int position);
	void removevertex(int tour, int position);
	void swapvertex(int tour, int i, int j);//assumption i < j
	void optvertices(int tour, int i, int j,bool breakreschedule);//assumption i < j

	friend ostream& operator<<(ostream& output, Sol& sol)
	{
		const char sep = ' ';
		for (int d = 0; d < (int) sol.tours.size(); ++d)
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
		for (int d = 0; d < (int) sol.tours.size(); ++d)
		{
			int end = (int) sol.tours[d].seq.size();
			output << left << setw(4) << setfill(sep) << "d" << d << "\n";
			output << left << setw(4) << setfill(sep) << "i";
			output << left << setw(4) << setfill(sep) << "vi";
			output << left << setw(7) << setfill(sep) << "sco";
			output << left << setw(7) << setfill(sep) << "vol";
			output << left << setw(7) << setfill(sep) << "wei";
			output << left << setw(8) << setfill(sep) << "dep";
			output << left << setw(8) << setfill(sep) << "shift"<< "\n";
			for (int i = 0; i < end; ++i)
			{
				output << left << setw(4) << setfill(sep) << i;
				output << left << setw(4) << setfill(sep) << sol.tours[d].seq[i]->index;
				output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->score;
				output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->volume;
				output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->weight;
				output << left << setw(8) << setfill(sep) << sol.tours[d].deptime[i];
				output << left << setw(8) << setfill(sep) << sol.tours[d].max_shift[i]<< "\n";
			}
		}
		output << "total score : " << sol.score;
		return output;
	}//end output operator
};