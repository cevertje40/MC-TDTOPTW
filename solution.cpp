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
		tours[t].max_shift.push_back(0.0);
		tours[t].wait_at.reserve(ins.maxvertices);
		tours[t].wait_at.push_back(0.0);
		tours[t].pref_wait.reserve(ins.maxvertices);
		tours[t].pref_wait.push_back(0.0);
		tours[t].br_margin.reserve(ins.maxvertices);
		tours[t].br_margin.push_back(0.0);
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

void Tour::update(int start, int end) // keep break fixed and update travel time and max_shift
{
	const double EDT = ins->t[index].EDT;
	const double B = ins->breakdur;
	const int end_depot_idx = ins->maxvertices - 1;

	double t = deptime[start] + EDT;

	for (int u = start; u < end - 1; ++u) 
	{
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];

		double arr = ins->arrival_time(o->con[p->index], t);

		if (action[u + 1]) 
		{ // break at p
			const bool endp = (p->index == end_depot_idx);
			if (!endp && arr < ins->breakstart) arr = ins->breakstart; // wait to breakstart (non-depot)
			arr += B;                                                 // take the break
		}

		if (arr < p->LTW[index]) arr = p->LTW[index]; // LTW after break if any
		arr += p->serv;

		// delta update of slack
		max_shift[u + 1] = (deptime[u + 1] + max_shift[u + 1]) - (arr - EDT);

		deptime[u + 1] = arr - EDT;
		t = arr;
	}
}

void Tour::update_break(int newbreakindex)
{
	const int n = (int)seq.size();
	if (n < 2) return;                         // nothing to do
	if (newbreakindex < 1 || newbreakindex > n - 1) return;  // invalid target

	// clear old mark (if any) and set the new one
	if (breakindex >= 1 && breakindex <= n - 1) action[breakindex] = 0;
	breakindex = newbreakindex;
	action[newbreakindex] = 1;

	// forward recompute with fixed break position (start-only semantics)
	const double EDT = ins->t[index].EDT;
	const double B = ins->breakdur;
	const int END_DEPOT = ins->maxvertices - 1;

	deptime[0] = 0.0;
	double t = EDT;

	for (int u = 0; u < n - 1; ++u)
	{
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];

		double arr = ins->arrival_time(o->con[p->index], t);

		// Take the break exactly at the chosen vertex (before serving p)
		if (u + 1 == breakindex)
		{
			const bool endp = (p->index == END_DEPOT);
			// Non-depot: if early, wait to breakstart; depot: allowed to start on arrival
			double start = endp ? arr : std::max(arr, ins->breakstart);
			arr = start + B;
		}

		// Enforce LTW after the break (if any), then add service
		if (arr < p->LTW[index]) arr = p->LTW[index];
		arr += p->serv;

		deptime[u + 1] = arr - EDT;
		t = arr;
	}

	// Recompute max_shift consistently with this fixed break position
	calc_maxshift();
}


bool Tour::update_break() {
	const int n = (int)seq.size();
	const double EDT = ins->t[index].EDT;
	const double LAT = ins->t[index].LAT;
	const double B = ins->breakdur;
	const double BEND = ins->breakend;
	const int   END = ins->maxvertices - 1;

	// 1) Build NO-BREAK forward locally
	std::vector<double> arr0(n, 0.0), svc(n, 0.0), dep_no_break(n, 0.0);
	double t = EDT;
	dep_no_break[0] = 0.0;
	for (int u = 0; u < n - 1; ++u) {
		auto* o = seq[u]; auto* p = seq[u + 1];
		double a = ins->arrival_time(o->con[p->index], t);
		arr0[u + 1] = a;
		double s = (a < p->LTW[index] ? p->LTW[index] : a);
		// If even the no-break start-of-service exceeds UTW, no feasible break will fix it.
		if (s > p->UTW[index] + 1e-9) return false;
		svc[u + 1] = s;
		t = s + p->serv;
		dep_no_break[u + 1] = t - EDT;
	}

	// 2) Find a feasible break start index
	auto can_place = [&](int u)->bool {
		double cur = EDT;  // recompute from start with the break at u
		for (int s = 0; s < n - 1; ++s) {
			auto* o = seq[s]; auto* p = seq[s + 1];
			double arr = ins->arrival_time(o->con[p->index], cur);
			if (s + 1 == u) {
				bool depot = (p->index == END);
				double latest = depot ? std::min(BEND, LAT - B)
					: std::min(BEND, p->UTW[index] - B);
				if (!depot && arr < ins->breakstart) arr = ins->breakstart;
				if (arr > latest + 1e-9) return false;
				arr += B;
			}
			if (arr < p->LTW[index]) arr = p->LTW[index];
			if (arr > p->UTW[index] + 1e-9) return false;
			cur = arr + p->serv;
		}
		return true;
		};

	int chosen = -1;
	for (int j = 1; j < n; ++j) {
		bool depot = (seq[j]->index == END);
		double earliest = depot ? arr0[j] : std::max(arr0[j], ins->breakstart);
		double latest = depot ? std::min(BEND, LAT - B)
			: std::min(BEND, seq[j]->UTW[index] - B);
		if (earliest <= latest + 1e-9) {
			if (can_place(j)) { chosen = j; break; }
		}
	}
	if (chosen < 0) return false;

	// 3) Commit atomically
	std::fill(action.begin(), action.end(), 0);
	action[chosen] = 1;
	breakindex = chosen;

	// 4) FULL forward recompute of deptime under break semantics
	double cur = EDT;
	deptime[0] = 0.0;
	for (int s = 0; s < n - 1; ++s) {
		auto* o = seq[s]; auto* p = seq[s + 1];
		double arr = ins->arrival_time(o->con[p->index], cur);
		if (s + 1 == chosen) {
			bool depot = (p->index == END);
			if (!depot && arr < ins->breakstart) arr = ins->breakstart;
			arr += B;
		}
		if (arr < p->LTW[index]) arr = p->LTW[index];
		cur = arr + p->serv;
		deptime[s + 1] = cur - EDT;
	}

	// 5) Refresh slack + budgets once
	calc_maxshift();
	return true;
}

void Tour::update_maxshift(int start, int end, double arrivaltime)
{
	// ===== 1) Backward cap for max_shift on [start+1 .. end]
	for (int i = end; i > start; --i)
	{
		Ins::Vertex* y = seq[i];
		Ins::Vertex* z = seq[i + 1];

		double L_dep_i = ins->departure_time(y->con[z->index], arrivaltime);

		// cap by y’s UTW + serv (latest start-of-service + service)
		if (L_dep_i > y->UTW[index] + y->serv) L_dep_i = y->UTW[index] + y->serv;

		// if break at y, departure (after break) ≤ breakend + breakdur
		if (action[i] == 1) {
			double cap = ins->breakend + ins->breakdur;
			if (L_dep_i > cap) L_dep_i = cap;
		}

		// slack at y (relative to scheduled departure EDT + deptime[i])
		max_shift[i] = L_dep_i - (ins->t[index].EDT + deptime[i]);

		// move the “latest arrival-before-service” cap to predecessor
		arrivaltime = L_dep_i - (y->serv + (action[i] ? ins->breakdur : 0.0));
	}

	// ===== 2) Refresh per-vertex budgets on a NO-BREAK forward pass
	const int n = (int)seq.size();
	const double EDT = ins->t[index].EDT;
	const double LAT = ins->t[index].LAT;
	const double B = ins->breakdur;
	const double BEND = ins->breakend;
	const int ENDDEP = ins->maxvertices - 1;

	wait_at.assign(n, 0.0);
	pref_wait.assign(n, 0.0);
	br_margin.assign(n, 0.0);

	// We also need transient arr0/svc to compute earliest and wait_at
	double t = EDT;
	std::vector<double> arr0(n, 0.0), svc(n, 0.0);

	for (int u = 0; u < n - 1; ++u)
	{
		auto* o = seq[u];
		auto* p = seq[u + 1];
		const bool depot = (p->index == ENDDEP);

		// arrival BEFORE any waiting (still NO break)
		double a = ins->arrival_time(o->con[p->index], t);
		arr0[u + 1] = a;

		// start-of-service under LTW only
		double s = (a < p->LTW[index] ? p->LTW[index] : a);
		svc[u + 1] = s;

		// usable waiting to START a break at this vertex
		double anchor = depot ? a : std::max(a, ins->breakstart);
		double w = s - anchor; if (w < 0.0) w = 0.0;
		wait_at[u + 1] = w;
		pref_wait[u + 1] = pref_wait[u] + w;

		// advance no-break clock
		t = s + p->serv;
	}

	// ===== 3) br_margin with downstream-slack cap (uses updated max_shift from step 1)
	for (int j = 1; j < n; ++j)
	{
		auto* p = seq[j];
		const bool depot = (p->index == ENDDEP);

		// earliest legal break start at boundary j (start-only semantics)
		double earliest = depot ? arr0[j] : std::max(arr0[j], ins->breakstart);

		// local latest bound (window & breakend)
		double latest_local = depot
			? std::min(BEND, LAT - B)
			: std::min(BEND, p->UTW[index] - B);

		// downstream-slack cap: we must still fit break + service before the latest allowed departure at j
		double dep_cap_j = EDT + deptime[j] + max_shift[j];       // latest feasible departure time *after* break/service
		double latest_from_slack = dep_cap_j - (p->serv + B);

		double latest = std::min(latest_local, latest_from_slack);
		double margin = latest - earliest;
		if (margin < 0.0) margin = 0.0;

		br_margin[j] = margin;
	}
}

void Tour::calc_maxshift() {
	const int n = (int)seq.size();
	const double EDT = ins->t[index].EDT;
	const double LAT = ins->t[index].LAT;
	const double TMAX = ins->t[index].T_max;
	const double B = ins->breakdur;
	const double BEND = ins->breakend;
	const int ENDDEP = ins->maxvertices - 1;

	wait_at.assign(n, 0.0);
	pref_wait.assign(n, 0.0);
	br_margin.assign(n, 0.0);

	// ---- 1) Backward cap on the CURRENT schedule (uses existing deptime/action)
	double L_dep_end = std::min(LAT, EDT + TMAX);
	if (action.back() == 1) L_dep_end = std::min(L_dep_end, BEND + B);
	max_shift.back() = L_dep_end - (EDT + deptime.back());

	double arr_cap = (action.back() == 1) ? (std::min(LAT, BEND + B) - B) : LAT;
	for (int i = n - 2; i > 0; --i) {
		auto* y = seq[i];
		auto* z = seq[i + 1];
		const bool has_break = (action[i] == 1);

		double L_dep_i = ins->departure_time(y->con[z->index], arr_cap);
		L_dep_i = std::min(L_dep_i, y->UTW[index] + y->serv);
		if (has_break) L_dep_i = std::min(L_dep_i, BEND + B);

		max_shift[i] = L_dep_i - (EDT + deptime[i]);
		arr_cap = L_dep_i - (y->serv + (has_break ? B : 0.0));
	}

	// ---- 2) Transient NO-BREAK forward pass for arr0/svc + wait_at/pref_wait
	std::vector<double> arr0(n, 0.0), svc(n, 0.0);
	double t = EDT;
	for (int u = 0; u < n - 1; ++u) {
		auto* o = seq[u];
		auto* p = seq[u + 1];

		const double a = ins->arrival_time(o->con[p->index], t);   // BEFORE waiting
		const double s = (a < p->LTW[index] ? p->LTW[index] : a);  // LTW only
		arr0[u + 1] = a;
		svc[u + 1] = s;

		// usable waiting at boundary u+1 (start-of-service − anchor)
		const bool depot = (p->index == ENDDEP);
		const double anchor = depot ? a : std::max(a, ins->breakstart);
		const double w = std::max(0.0, s - anchor);
		wait_at[u + 1] = w;
		pref_wait[u + 1] = pref_wait[u] + w;

		t = s + p->serv; // advance no-break clock
	}

	// ---- 3) Compute br_margin with a downstream-slack cap
	for (int j = 1; j < n; ++j) {
		auto* p = seq[j];
		const bool depot = (p->index == ENDDEP);

		const double earliest = depot ? arr0[j]
			: std::max(arr0[j], ins->breakstart);

		const double latest_local = depot
			? std::min(BEND, LAT - B)
			: std::min(BEND, p->UTW[index] - B);

		// latest departure (after break+service) allowed at boundary j
		const double dep_cap_j = EDT + deptime[j] + max_shift[j];

		// break start must leave time for break + service
		const double latest_from_slack = dep_cap_j - (p->serv + B);

		const double latest = std::min(latest_local, latest_from_slack);

		double margin = latest - earliest;
		if (margin < 0.0) margin = 0.0;
		br_margin[j] = margin;
	}
}

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
	update_break();
	calc_maxshift();
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

void Tour::opt_vertices(int i, int j,int newbreakindex)
{
	//reverse sequence
	for (int f = 0; f < 1 + (j - i) / 2; ++f)
	{
		Ins::Vertex* temp = seq[j - f];
		seq[j - f] = seq[i + f];
		seq[i + f] = temp;
	}
	if (newbreakindex != breakindex)
	{
		action[breakindex] = 0;
		action[newbreakindex] = 1;
		breakindex = newbreakindex;
		update(0, int(seq.size()));
		calc_maxshift();
	}
	else
	{
		update(i - 1, int(seq.size()));
		calc_maxshift();
	}
}

void Tour::swap_vertices(int i, int j,int newbreakindex)
{
	Ins::Vertex* remember = seq[i];
	seq[i] = seq[j];
	seq[j] = remember;
	if (newbreakindex!=breakindex)
	{
		action[breakindex] = 0;
		action[newbreakindex] = 1;
		breakindex = newbreakindex;
		update(0,int(seq.size()));
		calc_maxshift();
	}
	else
	{
		update(i - 1, int(seq.size()));
		calc_maxshift();
	}
}





bool Tour::check()
{
	bool tourok = true;
	int scorecheck = 0;
	double weightcheck = 0;
	double volumecheck = 0;

	const double EDT = ins->t[index].EDT;
	const double LAT = ins->t[index].LAT;
	const double B = ins->breakdur;
	const double BE = ins->breakend;
	const int ENDDEP = ins->maxvertices - 1;

	// a. travel time check
	double currenttime = EDT + deptime[0];
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

		double startbreak = 0.0;
		double endbreak = 0.0;
		double startservice = 0.0;
		double endservice = 0.0;

		if (breakcurrent)
		{
			const bool endp = (current->index == ENDDEP);

			// break can start upon arrival; at non-depot clamp to breakstart
			if (!endp && arrivaltime < ins->breakstart)
				startbreak = ins->breakstart;
			else
				startbreak = arrivaltime;

			// NEW: latest legal start guard
			double latest_start = endp
				? std::min(BE, LAT - B)
				: std::min(BE, current->UTW[index] - B);

			if (startbreak > latest_start + 1e-6)
			{
				std::cout << term::fg(term::Color::red)
					<< "tour: " << index
					<< " FAILURE!!! break start exceeds latest_start at pos=" << (i + 1)
					<< " node=" << current->index
					<< " start=" << startbreak << " latest=" << latest_start << "\n";
				tourok = false;
			}

			endbreak = startbreak + B;

			// Apply LTW at start-of-service
			startservice = (endbreak < current->LTW[index]) ? current->LTW[index] : endbreak;
		}
		else
		{
			// No break: just LTW
			startservice = (arrivaltime < current->LTW[index]) ? current->LTW[index] : arrivaltime;
		}

		// Windows at start-of-service
		if (startservice < current->LTW[index] - 1e-6)
		{
			std::cout << term::fg(term::Color::red)
				<< "tour: " << index << " FAILURE!!! LTW fail at pos=" << (i + 1)
				<< " node=" << current->index << "\n";
			tourok = false;
		}
		if (startservice > current->UTW[index] + 1e-6)
		{
			std::cout << term::fg(term::Color::red)
				<< "tour: " << index << " FAILURE!!! UTW fail at pos=" << (i + 1)
				<< " node=" << current->index << " start=" << startservice
				<< " UTW=" << current->UTW[index] << "\n";
			tourok = false;
		}

		endservice = startservice + current->serv;
		currenttime = endservice;

		// NEW: consistency check at each step
		double stored = EDT + deptime[i + 1];
		if (std::fabs(stored - currenttime) > 1e-6)
		{
			std::cout << term::fg(term::Color::red)
				<< "tour: " << index
				<< " deptime mismatch at pos=" << (i + 1)
				<< " calc=" << (currenttime - EDT)
				<< " stored=" << deptime[i + 1] << "\n";
			tourok = false;
		}
	}

	length = currenttime - EDT;

	// NEW: end depot LAT guard
	if (currenttime > LAT + 1e-6)
	{
		std::cout << term::fg(term::Color::red)
			<< "tour: " << index << " arrival " << currenttime
			<< " above LAT " << LAT << "\n";
		tourok = false;
	}

	if (length > ins->t[index].T_max + 1e-6)
	{
		std::cout << term::fg(term::Color::red)
			<< "tour: " << index << " length " << length
			<< " above T_max " << ins->t[index].T_max << "\n";
		tourok = false;
	}

	if (std::fabs(length - deptime.back()) > 1e-6)
	{
		std::cout << term::fg(term::Color::red)
			<< "tour: " << index << " total length mismatch. calc=" << length
			<< " stored=" << deptime.back() << " T_max=" << ins->t[index].T_max << "\n";
		tourok = false;
	}

	// c. weight checks
	if (std::fabs(weightcheck - weight) > 1e-6)
	{
		std::cout << term::fg(term::Color::red)	<< "tour: " << index << " new weight " << weightcheck<< " stored " << weight << " max " << ins->t[index].W_max << "\n";
		tourok = false;
	}
	if (weight - ins->t[index].W_max > 1e-6)
	{
		std::cout << term::fg(term::Color::red)	<< "tour: " << index << " weight " << weight<< " above max " << ins->t[index].W_max << "\n";
		tourok = false;
	}

	// d. volume checks
	if (std::fabs(volumecheck - volume) > 1e-6)
	{
		std::cout << term::fg(term::Color::yellow)<< "tour: " << index << " new volume " << volumecheck<< " stored " << volume << " max " << ins->t[index].V_max << "\n";
		tourok = false;
	}
	if (volume - ins->t[index].V_max > 1e-6)
	{
		std::cout << term::fg(term::Color::red)	<< "tour: " << index << " volume " << volume<< " above max " << ins->t[index].V_max << "\n";
		tourok = false;
	}

	// e. exactly one break, index consistency
	int amountbreaks = 0;
	for (int i = 0; i <= end; ++i) // break can be at end depot
	{
		if (action[i] == 1)
		{
			if (breakindex != i)
			{
				std::cout << term::fg(term::Color::red)
					<< "tour: " << index << " breakindex and action mismatch\n";
				tourok = false;
			}
			++amountbreaks;
		}
	}
	if (amountbreaks != 1)
	{
		std::cout << term::fg(term::Color::red)
			<< "tour: " << index << " amount of breaks not ok " << amountbreaks << "\n";
		tourok = false;
	}

	// g. start/end depots
	if (!((seq[0] == &ins->v[0]) && (seq.back() == &ins->v[ins->maxvertices - 1])))
	{
		std::cout << term::fg(term::Color::red)
			<< "tour: " << index << " start/end vertex not ok\n";
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
		tours[tour].max_shift.push_back(0.0);
		tours[tour].wait_at.clear();
		tours[tour].wait_at.push_back(0.0);
		tours[tour].pref_wait.clear();
		tours[tour].pref_wait.push_back(0.0);
		tours[tour].br_margin.clear();
		tours[tour].br_margin.push_back(0.0);
		tours[tour].score = 0;
		tours[tour].weight = 0;
		tours[tour].volume = 0;
		tours[tour].breakindex = -1;
	}
	score = 0;
	available.set();//sets all bits to true
	available[ins->v[0].index] = false;
}
