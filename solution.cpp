#include "solution.h"

using namespace std;

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

void Sol::read_from_file()
{
	vector<vector<int>> inputv(ins->maxtours);
	vector<vector<int>> inputb(ins->maxtours);
	ifstream ifs;
	ifs.open("debug_sol.txt", ifstream::in);
	if (ifs.is_open())
	{
		for (int t = 0; t < ins->maxtours; ++t)
		{
			int seqsize;
			string line;
			//read in number of vertices in tour
			getline(ifs, line);
			stringstream str(line);
			str >> seqsize;
			//read in vertex indices
			getline(ifs, line);
			str = stringstream(line);
			for (int i = 0; i < seqsize; ++i)
			{
				int index;
				str >> index;
				inputv[t].push_back(index);
			}
			//read in vertex action
			getline(ifs, line);
			str = stringstream(line);
			for (int i = 0; i < seqsize; ++i)
			{
				int action;
				str >> action;
				inputb[t].push_back(action);
			}
		}
		ifs.close();
	}
	else
	{
		printf("\ninput error in filenames file");
	}
	//construct the tours
	for (int t = 0; t < ins->maxtours; ++t)
	{
		Sol::Tour* tour = &tours[t];
		double currenttime = ins->t[t].EDT;
		Ins::Vertex* last = tour->seq[0];
		for (int i = 1; i < inputv[t].size(); ++i)
		{
			Ins::Vertex* candidate = &ins->v[inputv[t][i]];
			tour->seq.push_back(candidate);
			available[candidate->index] = false;
			tour->score += candidate->score;
			score += candidate->score;
			tour->max_shift.push_back(0);//dummy die dan in calc max shift upgedate wordt
			tour->volume += candidate->volume;
			tour->weight += candidate->weight;
			tour->action.push_back(inputb[t][i]);
			if (inputb[t][i] == 1)
			{
				tour->breakindex = i+1;
			}
			//deptime calc
			double arrivaltime = ins->arrival_time(last->con[candidate->index], currenttime);
			if (arrivaltime + inputb[t][i] * ins->breakdur < candidate->LTW[t])
			{
				arrivaltime = candidate->LTW[t] - (inputb[t][i] * ins->breakdur);
			}
			arrivaltime += candidate->serv + (inputb[t][i] * ins->breakdur);
			tour->deptime.push_back(arrivaltime - ins->t[t].EDT);
			currenttime = arrivaltime;
			last = candidate;
		}//end tour creation
		tour->calc_maxshift();
	}//end for all tours
}//end input custom

void Sol::write_to_file()
{
	ofstream output;
	output.open("debug_sol_test.txt", ios::out);
	for (int t = 0; t < ins->maxtours; ++t)
	{
		Sol::Tour* tour = &tours[t];
		output << tour->seq.size()<<"\n";
		for (int i = 0; i < tour->seq.size();++i)
		{
			output<<tour->seq[i]->index << " ";
		}
		output << "\n";
		for (int i = 0; i < tour->seq.size(); ++i)
		{
			output << tour->action[i] << " ";
		}
		output << "\n";
	}
	output.close();
}

void Sol::write_to_cplex()
{
	ofstream output;
	string name = "sol_" + to_string(ins->maxvertices)+ ".txt";
	output.open(name, ios::out);
	output << score << "\n";
	for (int t = 0; t < ins->maxtours; ++t)
	{
		Sol::Tour* tour = &tours[t];
		output << tour->seq.size()<<" ";
		for (int i = 0; i < tour->seq.size(); ++i)
		{
			output << tour->seq[i]->index << " ";
		}
		output << "\n";
		for (int i = 0; i < tour->seq.size(); ++i)
		{
			output << time_periods[0] + tour->deptime[i] << " ";
		}
		output << "\n";
		for (int i = 0; i < tour->seq.size(); ++i)
		{
			output << ins->find_t(time_periods[0] + tour->deptime[i])<<" ";
		}
		output << "\n";
		for (int i = 0; i < tour->seq.size(); ++i)
		{
			if (tour->action[i] == 1)
			{
				output << tour->seq[i]->index << endl;
			}
		}
		output << "\n";
	}
	output.close();
}

void Sol::check_availability()
{
	for (int i = 0;i < ins->maxvertices;++i)
	{
		if (available[i])
		{
			cout << "i: " << i << "still available" << endl;
		}
	}
}

int Sol::repair()
{
	int total_removed = 0;
	for (int t = 0; t < ins->maxtours; ++t) {
		Sol::Tour* tour = &tours[t];
		auto result = tour->repair();
		score -= result.first;      // subtract score decrease
		total_removed += result.second; // accumulate removed vertices
	}
	std::cout << "Vertices removed in repair: " << total_removed << std::endl;
	return total_removed;
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
		output << left << setw(4) << setfill(sep) << "sco";
		output << left << setw(7) << setfill(sep) << "vol";
		output << left << setw(7) << setfill(sep) << "wei";
		output << left << setw(8) << setfill(sep) << "dep";
		output << left << setw(8) << setfill(sep) << "LTW";
		output << left << setw(8) << setfill(sep) << "UTW";
		output << left << setw(8) << setfill(sep) << "shift";
		output << left << setw(4) << setfill(sep) << "break" << "\n";
		for (int i = 0; i < end; ++i)
		{
			output << left << setw(4) << setfill(sep) << i;
			output << left << setw(4) << setfill(sep) << sol.tours[d].seq[i]->index;
			output << left << setw(4) << setfill(sep) << sol.tours[d].seq[i]->score;
			output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->volume;
			output << left << setw(7) << setfill(sep) << sol.tours[d].seq[i]->weight;
			output << left << setw(8) << setfill(sep) << sol.tours[d].deptime[i];
			output << left << setw(8) << setfill(sep) << sol.tours[d].seq[i]->LTW[d] - sol.ins->t[d].EDT;
			output << left << setw(8) << setfill(sep) << sol.tours[d].seq[i]->UTW[d] - sol.ins->t[d].EDT;
			output << left << setw(8) << setfill(sep) << sol.tours[d].max_shift[i];
			output << left << setw(4) << setfill(sep) << sol.tours[d].action[i] << "\n";
		}
	}
	output << "total score : " << sol.score;
	return output;
}//end output operator

bool Sol::check()
{
	//1. check tour constraints
	bool solok = true;
	int scorecheck = 0;
	vector<int> included(ins->maxvertices, 0);
	for (int d = 0; d < ins->maxtours; ++d)
	{
		if (!tours[d].check())
		{
			solok = false;
		}
		for (int i = 0; i < tours[d].seq.size(); ++i)
		{
			++included[tours[d].seq[i]->index];
		}
		scorecheck += tours[d].score;
	}//end for all tours
	//2. inclusion & availability check
	for (int i = 1; i < ins->maxvertices-1; ++i)//depot can be included multiple times
	{
		if (included[i] > 1)//regular vertex included more than once
		{
			cout << term::fg(term::Color::red) << "included more than one, vertex:" << i << endl;
			solok = false;
		}
		if (included[i] == 1)
		{
			if (available[i] == true)
			{
				cout << term::fg(term::Color::red) << "error in availability bitset for vertex: " << i << endl;
				solok = false;
			}
		}
		if (included[i] == 0)
		{
			if (available[i] == false)
			{
				cout << term::fg(term::Color::red) << "error in availability bitset for vertex: " << i << endl;
				solok = false;
			}
		}
	}
	//3. total score check
	if (scorecheck != score)
	{
		cout << term::fg(term::Color::red) << "score of best solution should be: " << scorecheck << " stored score is: " << score << endl;
		solok = false;
	}
	return solok;
}

void Tour::update(int start, int end)//keep break fixed and update travel time and max_shift
{
	//update after start
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

void Tour::update_break(int newbreakindex)
{
	//update break to new position
	action[breakindex] = 0;
	breakindex = newbreakindex;
	action[newbreakindex] = 1;
	//update complete solution because break is repositioned
	double currenttime = ins->t[index].EDT;
	int end = (int)seq.size();
	for (int u = 0; u < end - 1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if (arrivaltime + (action[u + 1] * ins->breakdur) < p->LTW[index])
		{
			arrivaltime = p->LTW[index] - (action[u + 1] * ins->breakdur);
		}
		arrivaltime += p->serv + (action[u + 1] * ins->breakdur);
		deptime[u + 1] = arrivaltime - ins->t[index].EDT;
		currenttime = arrivaltime;
	}//end for all u
	calc_maxshift();
}

void Tour::update_break()  // reposition break, update times and max_shift
{
	/*
	const int n = (int)seq.size();
	// reset actions
	std::fill(action.begin(), action.end(), 0);
	breakindex = -1;

	// --- forward pass WITHOUT break: record arr0 (arrival before any wait), svc (start-of-service) ---
	std::vector<double> arr0(n, 0.0);
	std::vector<double> svc(n, 0.0);

	double t = ins->t[index].EDT;
	deptime[0] = 0.0;

	for (int u = 0; u < n - 1; ++u) {
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];
		double a = ins->arrival_time(o->con[p->index], t); // arrival BEFORE waiting
		arr0[u + 1] = a;

		double s = (a < p->LTW[index] ? p->LTW[index] : a); // apply LTW wait
		svc[u + 1] = s;

		double leave = s + p->serv;   // (no break in this phase)
		deptime[u + 1] = leave - ins->t[index].EDT;
		t = leave;
	}

	// max_shift on the no-break schedule
	calc_maxshift();

	// --- usable waiting AFTER breakstart at each vertex (this is the key fix) ---
	//     wait_after_start[j] = max(0, svc[j] - max(arr0[j], breakstart))
	std::vector<double> wait_after_start(n, 0.0);
	for (int j = 1; j < n; ++j) {
		double usable = svc[j] - std::max(arr0[j], ins->breakstart);
		if (usable < 0.0) usable = 0.0;
		wait_after_start[j] = usable;
	}

	// --- earliest feasible position to start a break (start-only semantics) ---
	//     use your original idea, but with arr0 instead of current-time t
	int earliestbreakindex = n - 1; // default to end depot
	for (int u = 0; u < n - 1; ++u) {
		// We can start a break at p=seq[u+1] if either:
		// - we arrive after breakstart, OR
		// - we could wait-to-LTW and start break at LTW - breakdur (>= breakstart).
		//   That is: max(arr0[u+1], LTW[p] - B) >= breakstart
		double lhs = std::max(arr0[u + 1], seq[u + 1]->LTW[index] - ins->breakdur);
		if (lhs >= ins->breakstart - 1e-9) {
			earliestbreakindex = u + 1;
			break;
		}
	}

	// --- choose the first position that can absorb the break ---
	bool placed = false;
	for (int u = earliestbreakindex; u <= n - 1; ++u) {     // end depot included
		// IMPORTANT: use usable waiting after breakstart, not total early wait
		if (max_shift[u] + wait_after_start[u] + 1e-9 >= ins->breakdur) {
			action[u] = 1;
			breakindex = u;
			placed = true;
			break;
		}
	}
	if (!placed) {
		// fallback: put break on end depot (your original special case)
		action.back() = 1;
		breakindex = n - 1;
	}

	// --- recompute forward from breakindex-1 with the inserted break flags ---
	double currenttime = (breakindex > 0)
		? deptime[breakindex - 1] + ins->t[index].EDT
		: ins->t[index].EDT;

	for (int u = std::max(0, breakindex - 1); u < n - 1; ++u) {
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];

		double arr = ins->arrival_time(o->con[p->index], currenttime);

		// If there is a break at p, you can start it upon arrival, but if p is not the end depot
		// and you arrive before breakstart, the extra waiting-to-breakstart happens here.
		if (action[u + 1]) {
			if (p->index != ins->maxvertices - 1 && arr < ins->breakstart)
				arr = ins->breakstart;                    // eat only "usable" wait
			arr += ins->breakdur;                         // break finishes
		}

		if (arr < p->LTW[index]) arr = p->LTW[index];     // LTW wait after break if needed

		arr += p->serv;
		// update max_shift relative to the old time
		max_shift[u + 1] = (deptime[u + 1] + max_shift[u + 1]) - (arr - ins->t[index].EDT);
		if (max_shift[u + 1] < 0.0) max_shift[u + 1] = 0.0;

		deptime[u + 1] = arr - ins->t[index].EDT;
		currenttime = arr;
	}

	// --- tiny safety: recompute end-depot time from its predecessor ---
	{
		double tp = ins->t[index].EDT + deptime[n - 2];
		double arr = ins->arrival_time(seq[n - 2]->con[seq[n - 1]->index], tp);
		if (action[n - 1]) arr += ins->breakdur;          // if break at depot
		deptime[n - 1] = arr - ins->t[index].EDT;
	}

	// --- end-depot slack and backward max_shift cap (unchanged from your code) ---
	if (breakindex == n - 1) 
	{
		// break at end depot
		max_shift.back() = (ins->t[index].T_max - deptime.back());
		double arrive_cap = ins->t[index].LAT - ins->breakdur;
		if (ins->t[index].EDT + ins->t[index].T_max > ins->breakend + ins->breakdur) {
			max_shift.back() = (ins->breakend + ins->breakdur) - (deptime.back() + ins->t[index].EDT);
			arrive_cap = ins->breakend;
		}
		update_maxshift(0, n - 2, arrive_cap);
	}
	else 
	{
		// break before depot → cap by the successor’s latest arrival-before-service
		double cap_at_next = (ins->t[index].EDT + deptime[breakindex + 1] + max_shift[breakindex + 1])
			- (seq[breakindex + 1]->serv + action[breakindex + 1] * ins->breakdur);
		update_maxshift(0, breakindex, cap_at_next);
	}
	*/
	//evaluate tour without break
	int end = (int)seq.size();
	for (int vv = 0; vv < action.size(); ++vv)
	{
		action[vv] = 0;
	}
	double currenttime = ins->t[index].EDT;
	int earliestbreakindex = (int)seq.size() - 1;//default to end depot
	bool earliestbreaktofind = true;
	vector<double>waitingtime(seq.size(), 0.0);
	for (int u = 0; u < end - 1; ++u)
	{
		//gather departure time and corresponding time slot
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];
		//travel time from van o to p
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		if ((earliestbreaktofind) && (max(p->LTW[index] - ins->breakdur, arrivaltime) >= ins->breakstart))//arrival time need to be higher than breakstart not dep-serv
		{//determine earliest feasible point to insert a break
			earliestbreakindex = u + 1;
			earliestbreaktofind = false;
		}
		if (arrivaltime < p->LTW[index])
		{
			waitingtime[u + 1] = p->LTW[index] - arrivaltime;
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
	for (int u = earliestbreakindex; u <= end - 1; ++u)//end depot included
	{
		if (max_shift[u] + waitingtime[u] >= ins->breakdur)
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
		breakindex = (int)seq.size() - 1;
		//cout << "forced break on end depot" << endl;

	}
	if (breakindex == seq.size() - 1)
	{
		//update redeparture time at end depot
		currenttime = deptime[breakindex - 1] + ins->t[index].EDT;
		Ins::Vertex* o = seq[breakindex - 1];
		Ins::Vertex* p = seq[breakindex];
		double arrivaltime = ins->arrival_time(o->con[p->index], currenttime);
		arrivaltime += ins->breakdur;
		deptime.back() = arrivaltime - ins->t[index].EDT;
		//update max_shift of end depot and all preceding vertices
		max_shift.back() = (ins->t[index].T_max - deptime.back());
		arrivaltime = ins->t[index].LAT - ins->breakdur;
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
	// Initialize max_shift for the end depot
	int size = (int)seq.size();
	double departuretime = 0.0;
	max_shift.back() = (ins->t[index].T_max - deptime.back());
	// Calculate arrival time, adjusting for break duration if necessary
	double arrivaltime = ins->t[index].LAT - (action.back() * ins->breakdur);//if you break at the end depot subtract breakduration
	// Adjust max_shift and arrivalTime if there's a break at the end depot
	if ((action.back() == 1)&&(ins->t[index].LAT > ins->breakend + ins->breakdur))
	{
		max_shift.back() = (ins->breakend + ins->breakdur) - (deptime.back() + ins->t[index].EDT);
		arrivaltime = ins->breakend;//remove breakduration
	}
	// Iterate backward through the sequence (excluding depots)
	for (int i = size-2; i >0; --i)
	{
		// Calculate the departure time
		Ins::Vertex* y = seq[i];
		int breaki = action[i];
		Ins::Vertex* z = seq[i+1];
		departuretime = ins->departure_time(y->con[z->index], arrivaltime);
		// Time window check on departure time
		if (departuretime > y->UTW[index] + y->serv)
		{
			departuretime = y->UTW[index] + y->serv;
		}
		// Adjust departure time if a break is scheduled
		if ((breaki == 1)&&(departuretime > ins->breakend + ins->breakdur))
		{
			departuretime = ins->breakend + ins->breakdur;
		}
		// Store the maximum shift result
		max_shift[i] = departuretime - (deptime[i] + ins->t[index].EDT);
		// Update arrivalTime for the next iteration
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
	update(position, int(seq.size()));//update travel time and maxshift for all positions after insertion
	double arrivaltime = (deptime[position + 2] + ins->t[index].EDT + max_shift[position + 2]) - (seq[position + 2]->serv + action[position + 2] * ins->breakdur);//service time eraftrekken
	update_maxshift(0, position + 1, arrivaltime);//update maxshift for all positions before insertion
}

void Tour::remove_vertex(int position)
{
	bool reqbreak = false;
	if (position <= breakindex)
	{
		reqbreak = true;
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
		update_break();
	}
	else
	{
		update(position - 1, int(seq.size()));//update travel time for all vertices after deletion
		calc_maxshift();
	}
}//end remove_vertex tour version

void Tour::remove_vertices(int position1, int position2)
{
	if (position1 == position2) return;         // or assert
	if (position1 > position2) swap(position1, position2);
	Ins::Vertex* candidate1 = seq[position1];
	Ins::Vertex* candidate2 = seq[position2];
	score -= candidate1->score;// update score of the tour
	score -= candidate2->score;// update score of the tour
	volume -= candidate1->volume;
	weight -= candidate1->weight;
	volume -= candidate2->volume;
	weight -= candidate2->weight;
	//remove position 2 first (higher index)
	seq.erase(seq.begin() + position2);//insert point y after x
	deptime.erase(deptime.begin() + position2);//insert temporary value
	max_shift.erase(max_shift.begin() + position2);
	action.erase(action.begin() + position2);//insert regular visit action change later when necessary
	seq.erase(seq.begin() + position1);//insert point y after x
	deptime.erase(deptime.begin() + position1);//insert temporary value
	max_shift.erase(max_shift.begin() + position1);
	action.erase(action.begin() + position1);//insert regular visit action change later when necessary
	update_break();
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
		update_break();//update travel time and max_shift for complete tour
	}
	else
	{
		update(position - 1, int(seq.size()));//update travel time and maxshift for all positions after replacement
		double arrivaltime = (deptime[position + 1] + ins->t[index].EDT + max_shift[position + 1]) - (seq[position + 1]->serv + action[position + 1] * ins->breakdur);//service time eraftrekken
		update_maxshift(0, position, arrivaltime);//update maxshift for all positions before replacement
	}
}

void Tour::replace_vertex(Ins::Vertex* candidate, int position, int breakindex)
{
	Ins::Vertex* old = seq[position];
	seq[position] = candidate;//replace point old with candidate
	max_shift[position] = 0;//dummy value
	volume -= old->volume - candidate->volume;
	weight -= old->weight - candidate->weight;
	score += candidate->score - old->score;// update score of the new solution
	update_break(breakindex);//update travel time and max_shift for complete tour
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
		update_break();
	}
	else
	{
		update(i - 1, int(seq.size()));
		calc_maxshift();
	}
	
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
		update_break();
	}
	else
	{
		update(i - 1, int(seq.size()));
		calc_maxshift();
	}
}

vector<double> Tour::compute_wait_suffix()
{
	const int n = (int)seq.size();
	std::vector<double> wait_at(n, 0.0);      // actual waiting to satisfy LTW
	std::vector<double> usable_at(n, 0.0);    // part of waiting that can absorb ripple
	std::vector<double> wait_suffix(n, 0.0);

	const double EDT = ins->t[index].EDT;
	const double B = ins->breakdur;
	const double BSTART = ins->breakstart;

	// depart from start; service at seq[0] already included in deptime[0]
	double t = deptime[0] + EDT;

	for (int j = 0; j < n - 1; ++j)
	{
		Ins::Vertex* x = seq[j];
		Ins::Vertex* z = seq[j + 1];
		const int breakz = action[j + 1];

		// raw arrival before any waiting or break addition
		double arr_raw = ins->arrival_time(x->con[z->index], t);

		// actual schedule: enforce LTW (accounting for break if it is at z)
		double arr = arr_raw;
		double w = 0.0;
		if (arr + breakz * B < z->LTW[index]) {
			w = z->LTW[index] - (arr + breakz * B);
			arr = z->LTW[index] - breakz * B;
		}
		wait_at[j + 1] = w;

		// usable waiting: if a break is here, the portion before breakstart can’t absorb ripple
		if (breakz == 1) {
			// time we’d need just to reach breakstart (if we arrived too early)
			double unusable = std::max(0.0, BSTART - arr_raw);
			usable_at[j + 1] = std::max(0.0, w - unusable);
		}
		else {
			usable_at[j + 1] = w;
		}

		// advance to next vertex time (service + possible break at z)
		t = arr + z->serv + breakz * B;
	}

	// suffix sum of usable waiting
	double acc = 0.0;
	for (int i = n - 1; i >= 0; --i) 
	{
		acc += usable_at[i];
		wait_suffix[i] = acc;
	}
	return wait_suffix;
}

bool Tour::check()
{
	bool tourok = true;
	int scorecheck = 0;
	double weightcheck = 0;
	double volumecheck = 0;
	//a. travel time check
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
		double waitingtime = 0.0;
		double startbreak = 0.0;
		double endbreak = 0.0;
		double startservice = 0.0;
		double endservice = 0.0;
		if (breakcurrent) 
		{
			if ((arrivaltime < ins->breakstart)&&(current->index != ins->maxvertices - 1))
			{//break can not start before breakstart
				waitingtime += ins->breakstart - arrivaltime;
				startbreak = ins->breakstart;
			}
			else
			{//break can start upon arrival or before breakstart at enddepot
				startbreak = arrivaltime;
			}
			if (startbreak > ins->breakend + 1e-6)
			{
				cout << "tour: " << index << "FAILURE!!! break end time exceed breakend for solutionnr: " << i + 1 << " /vertex index: " << current->index << endl;
				tourok = false;
			}
			endbreak = startbreak + ins->breakdur;
			if (endbreak < current->LTW[index])
			{
				waitingtime += current->LTW[index] - endbreak;
				//cout << "waiting time for: "<<"i"<<i+1<<" , " <<seq[i + 1]->index << " <=> " << waitingtime << endl;
				startservice = current->LTW[index];
			}
			else
			{
				startservice = endbreak;
			}
		}
		else
		{
			if (arrivaltime < current->LTW[index])
			{
				waitingtime += current->LTW[index] - arrivaltime;
				//cout << "waiting time for: "<<"i"<<i+1<<" , " <<seq[i + 1]->index << " <=> " << waitingtime << endl;
				startservice = current->LTW[index];
			}
			else
			{
				startservice = arrivaltime;
			}
		}
		
		if (startservice < current->LTW[index]- 1e-6)
		{
			cout << term::fg(term::Color::red) << "tour: " << index << "FAILURE!!! LTW fail for solutionnr: " << i + 1 << " /vertex index: " << current->index << endl;
			tourok = false;
		}

		if (startservice > current->UTW[index]+ 1e-6)
		{
			cout << term::fg(term::Color::red) << "tour: " << index << "FAILURE!!! UTW fail for solutionnr: " << i + 1 << " /vertex index: " << current->index << endl;
			tourok = false;
		}
		endservice = startservice + current->serv;
		//cout<<i+1<<" calc traveltime: " << endservice-ins->t[d].EDT << " stored: " << tours[d].deptime[i + 1] << endl;
		currenttime = endservice;
	}//end for i
	length = currenttime - ins->t[index].EDT;
	if (length > ins->t[index].T_max + 1e-6)
	{
		cout << term::fg(term::Color::red) << "tour: " << index << "length " << length << " above max: " << ins->t[index].T_max << endl;
		tourok = false;

	}
	if (fabs(length - deptime.back()) > 1e-6)
	{
		cout << term::fg(term::Color::red) << "tour: " << index << "new calculated length: " << length << " stored length: " << deptime.back() << "max length" << ins->t[index].T_max << endl;
		tourok = false;
	}
	//c. weight checks
	if (fabs(weightcheck - weight) > 1e-6)
	{
		cout << term::fg(term::Color::red) << "tour: " << index << "new calculated weight" << weightcheck << "stored weight: " << weight << "max: " << ins->t[index].W_max << endl;
		tourok = false;
	}
	if (weight-ins->t[index].W_max> 1e-6)
	{
		tourok = false;
		cout << term::fg(term::Color::red) << "tour: " << index << "weight " << weight << " above max: " << ins->t[index].W_max << endl;
	}
	//d. volume checks
	if (fabs(volumecheck - volume) > 1e-6)
	{
		cout << term::fg(term::Color::yellow) << "tour: " << index << "new calculated volume" << volumecheck << "stored volume " << volume << "max: " << ins->t[index].V_max << endl;
		tourok = false;
	}
	if (volume - ins->t[index].V_max> 1e-6)
	{
		tourok = false;
		cout << term::fg(term::Color::red) << "tour: " << index << " volume " << volume << " above max: " << ins->t[index].V_max << endl;
	}
	//e. break timing check
	bool breakcheck = false;
	int amountbreaks = 0;
	for (int i = 0; i <= end; ++i)//break kan op enddepot zitten
	{
		if (action[i] == 1)
		{
			if (breakindex != i)
			{
				cout << term::fg(term::Color::red) << "tour: " << index << "breakindex and sol action don't match" << endl;
				tourok = false;
			}
			++amountbreaks;
		}
	}
	//f. amount of breaks check
	if (amountbreaks == 1)
	{
		//cout << yellow << "break ok" << endl;
	}
	else
	{
		cout << term::fg(term::Color::red) << "tour: " << index << "amount of breaks not ok" << amountbreaks << endl;
		tourok = false;
	}
	//g. tour should start and end at the respective depots
	if ((seq[0] == &ins->v[0]) && (seq.back() == &ins->v[ins->maxvertices - 1]))
	{
		//cout<<yellow << "start and end ok" << endl;
	}
	else
	{
		cout << term::fg(term::Color::red)  <<"tour: " << index << "start and end vertex not ok" << endl;
		tourok = false;
	}
	//h. max_shift check
	vector<double> max_shiftcheck(max_shift.size(), 0);
	int size = (int)seq.size();
	double departuretime = 0;
	max_shiftcheck.back() = (ins->t[index].T_max - deptime.back());
	double arrivaltime = ins->t[index].LAT - (action.back() * ins->breakdur);//if you break at the end depot subtract breakduration
	if (action.back() == 1)//break op enddepot
	{
		if (ins->t[index].LAT > ins->breakend + ins->breakdur)
		{
			//cout<<"path: "<<d<< " bij calc maxshift break op enddepot verhindert een maxshift: " << endl;
			max_shiftcheck.back() = (ins->breakend + ins->breakdur) - (ins->t[index].EDT+deptime.back());
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
		if (fabs(max_shift[i] - max_shiftcheck[i]) > 0.01)
		{
			cout << term::fg(term::Color::red) <<"tour: "<< index << "error in max_shift for position: " << i <<" new max_shift: "<<max_shiftcheck[i]<<" stored max_shift: "<<max_shift[i] << endl;
			tourok = false;
		}
	}
	return tourok;
}//end tour check

pair<int, int> Tour::repair()
{
	int scoredecrease = 0;
	int removed = 0;

	bool infeasible = false;
	double currenttime = ins->t[index].EDT + deptime[0];
	int end = static_cast<int>(seq.size()) - 1;

	// first feasibility check
	for (int i = 0; i < end; ++i) {
		Ins::Vertex* last = seq[i];
		Ins::Vertex* current = seq[i + 1];
		int breakcurrent = action[i + 1];
		double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);

		if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[index]) {
			arrivaltime = current->LTW[index] - breakcurrent * ins->breakdur;
		}
		if (arrivaltime > current->UTW[index]) {
			infeasible = true;
			break;
		}
		arrivaltime += current->serv + breakcurrent * ins->breakdur;
		currenttime = arrivaltime;
	}
	if (currenttime > ins->t[index].EDT + ins->t[index].T_max) {
		infeasible = true;
	}

	// repair loop
	while (infeasible) {
		end = static_cast<int>(seq.size()) - 1;
		if (end <= 1) break; // no regular vertices left

		scoredecrease += seq[end - 1]->score;
		++removed;
		remove_vertex(end - 1);

		// recheck feasibility
		infeasible = false;
		currenttime = ins->t[index].EDT + deptime[0];
		end = static_cast<int>(seq.size()) - 1;

		for (int i = 0; i < end; ++i) {
			Ins::Vertex* last = seq[i];
			Ins::Vertex* current = seq[i + 1];
			int breakcurrent = action[i + 1];
			double arrivaltime = ins->arrival_time(last->con[current->index], currenttime);

			if (arrivaltime + breakcurrent * ins->breakdur < current->LTW[index]) {
				arrivaltime = current->LTW[index] - breakcurrent * ins->breakdur;
			}
			if (arrivaltime > current->UTW[index]) {
				infeasible = true;
				break;
			}
			arrivaltime += current->serv + breakcurrent * ins->breakdur;
			currenttime = arrivaltime;
		}
		if (currenttime > ins->t[index].EDT + ins->t[index].T_max) {
			infeasible = true;
		}
	}

	return { scoredecrease, removed };
}

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

void Sol::replace_vertex(Tour& tour, Ins::Vertex* candidate, int position, int breakindex)
{
	Ins::Vertex* old = tour.seq[position];
	available[candidate->index] = false;
	available[old->index] = true;
	score += candidate->score - old->score;// update score of the new solution
	tour.replace_vertex(candidate, position, breakindex);
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
