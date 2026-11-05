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

void Sol::read_from_file(const std::string& dataset)
{
	struct RouteData
	{
		std::vector<int> vertices;
		std::vector<double> departures;
		std::vector<int> slots;
		int breakvertex = -1;
	};

	const std::string filename = "sol_" + dataset;
	std::ifstream ifs(filename);

	if (!ifs.is_open())
	{
		std::cerr << "Could not open solution file: " << filename << '\n';
		return;
	}

	double fileObjective = 0.0;
	if (!(ifs >> fileObjective))
		throw std::runtime_error("Could not read objective from " + filename);

	std::vector<RouteData> routes(ins->maxtours);

	// Parse the whole file before modifying the current solution.
	for (int t = 0; t < ins->maxtours; ++t)
	{
		int seqsize = 0;
		if (!(ifs >> seqsize) || seqsize < 2)
			throw std::runtime_error("Invalid route size in " + filename);

		RouteData& route = routes[t];
		route.vertices.resize(seqsize);
		route.departures.resize(seqsize - 1);
		route.slots.resize(seqsize - 1);

		for (int& vertex : route.vertices)
		{
			if (!(ifs >> vertex))
				throw std::runtime_error("Could not read route vertices from " + filename);
		}

		for (double& departure : route.departures)
		{
			if (!(ifs >> departure))
				throw std::runtime_error("Could not read departure times from " + filename);
		}

		for (int& slot : route.slots)
		{
			if (!(ifs >> slot))
				throw std::runtime_error("Could not read time-slot indices from " + filename);
		}

		if (!(ifs >> route.breakvertex))
			throw std::runtime_error("Could not read break vertex from " + filename);
	}

	const int endDepot = ins->maxvertices - 1;

	// Validate depot endpoints, break locations, and customer uniqueness.
	std::vector<bool> seenCustomer(ins->maxvertices, false);

	for (int t = 0; t < ins->maxtours; ++t)
	{
		const RouteData& route = routes[t];

		if (route.vertices.front() != 0 || route.vertices.back() != endDepot)
			throw std::runtime_error("A route does not start and end at the depots");

		int breakPositions = 0;
		for (int pos = 0; pos < static_cast<int>(route.vertices.size()); ++pos)
		{
			const int vertex = route.vertices[pos];

			if (vertex < 0 || vertex >= ins->maxvertices)
				throw std::runtime_error("Invalid vertex index in " + filename);

			if (vertex == route.breakvertex)
				++breakPositions;

			if (pos > 0 && pos < static_cast<int>(route.vertices.size()) - 1)
			{
				if (seenCustomer[vertex])
					throw std::runtime_error("A customer occurs more than once in " + filename);

				seenCustomer[vertex] = true;
			}
		}

		if (breakPositions != 1)
			throw std::runtime_error("Break vertex must occur exactly once on its route");
	}

	// Reset solution and availability state.
	score = 0;
	std::fill(available.begin(), available.end(), true);

	for (int t = 0; t < ins->maxtours; ++t)
	{
		Tour& tour = tours[t];
		const RouteData& route = routes[t];
		const int n = static_cast<int>(route.vertices.size());
		const double EDT = ins->t[t].EDT;

		tour.seq.clear();
		tour.action.assign(n, 0);
		tour.deptime.clear();
		tour.max_shift.assign(n, 0.0);

		tour.score = 0;
		tour.weight = 0.0;
		tour.volume = 0.0;
		tour.breakindex = -1;

		// Build the route and mark its break location.
		for (int pos = 0; pos < n; ++pos)
		{
			const int vertex = route.vertices[pos];
			tour.seq.push_back(&ins->v[vertex]);

			if (vertex == route.breakvertex)
			{
				tour.action[pos] = 1;
				tour.breakindex = pos;
			}

			if (pos > 0 && pos < n - 1)
			{
				available[vertex] = false;
				tour.score += ins->v[vertex].score;
				tour.weight += ins->v[vertex].weight;
				tour.volume += ins->v[vertex].volume;
				score += ins->v[vertex].score;
			}
		}

		// Replay the route using the same timing rules as Tour::check().
		tour.deptime.push_back(0.0);
		double currenttime = EDT;

		for (int pos = 1; pos < n; ++pos)
		{
			Ins::Vertex* previous = tour.seq[pos - 1];
			Ins::Vertex* current = tour.seq[pos];

			const double arrival = ins->arrival_time(
				previous->con[current->index],
				currenttime);

			double startservice = 0.0;

			if (tour.action[pos] == 1)
			{
				const bool atEndDepot = (current->index == endDepot);

				// At customers, wait until brk_start if arriving earlier.
				// At the end depot, the break may start upon arrival.
				double startbreak = arrival;
				if (!atEndDepot && startbreak < ins->breakstart)
					startbreak = ins->breakstart;

				if (startbreak > ins->breakend + 1e-6)
					throw std::runtime_error(
						"Break starts after brk_end on route " + std::to_string(t));

				const double endbreak = startbreak + ins->breakdur;
				startservice = std::max(endbreak, current->LTW[t]);
			}
			else
			{
				startservice = std::max(arrival, current->LTW[t]);
			}

			if (startservice > current->UTW[t] + 1e-6)
				throw std::runtime_error(
					"Latest service time violated on route " + std::to_string(t));

			currenttime = startservice + current->serv;
			tour.deptime.push_back(currenttime - EDT);
		}

		tour.calc_maxshift();
	}

}

void Sol::write_to_file(const std::string& dataset)
{
	std::ofstream output;
	std::string name = "sol_" + dataset;   
	output.open(name, std::ios::out);

	if (!output.is_open())
	{
		std::cerr << "Could not open " << name << " for writing.\n";
		return;
	}

	output << score << "\n";

	for (int t = 0; t < ins->maxtours; ++t)
	{
		Sol::Tour* tour = &tours[t];
		const int size = static_cast<int>(tour->seq.size());

		// 1) route size
		output << size << " ";

		// 2) vertex sequence
		for (int i = 0; i < size; ++i)
		{
			output << tour->seq[i]->index << " ";
		}
		output << "\n";

		// 3) departure times for arcs: one value per arc (size - 1 values)
		// deptime[a] is the departure/service completion time at seq[a],
		// corresponding to arc (seq[a], seq[a+1]).
		for (int a = 0; a < size - 1; ++a)
		{
			output << tour->deptime[a] << " ";
		}
		output << "\n";

		// 4) time slots for arcs: one value per arc (size - 1 values)
		for (int a = 0; a < size - 1; ++a)
		{
			double abs_departure = ins->t[t].EDT + tour->deptime[a];
			output << ins->find_t(abs_departure) << " ";
		}
		output << "\n";

		// 5) one break vertex index
		// CPLEX expects q[p][breakvertex[p]] = 1
		int breakvertex = 0;
		if (tour->breakindex >= 0 && tour->breakindex < size)
		{
			breakvertex = tour->seq[tour->breakindex]->index;
		}
		output << breakvertex << "\n";
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


SolutionPerfStats Sol::repair_and_collect_stats()
{
	SolutionPerfStats stats;
	stats.route_stats.reserve(ins->maxtours);

	int removed_score_sum = 0;
	double removed_tw_width_sum = 0.0;
	double removed_service_sum = 0.0;
	double removed_weight_sum = 0.0;
	double removed_volume_sum = 0.0;
	double removed_depot_tt_sum = 0.0;
	double removed_position_sum = 0.0;

	double sum_route_duration = 0.0;
	double sum_route_time_util = 0.0;
	double sum_weight_util = 0.0;
	double sum_volume_util = 0.0;
	double sum_wait_per_customer = 0.0;
	double sum_break_start = 0.0;
	double sum_break_pos = 0.0;

	int routes_with_break_start = 0;
	int routes_with_break_pos = 0;
	int routes_break_at_end = 0;

	for (int t = 0; t < ins->maxtours; ++t)
	{
		Tour* tour = &tours[t];

		RemovedVertexStats removed_stats;
		auto result = tour->repair(removed_stats);

		score -= result.first;

		RoutePerfStats r = tour->collect_route_stats();
		r.removed = removed_stats;
		stats.route_stats.push_back(r);

		stats.total_removed_customers += removed_stats.removed_count;

		removed_score_sum += removed_stats.removed_score_sum;
		removed_tw_width_sum += removed_stats.removed_tw_width_sum;
		removed_service_sum += removed_stats.removed_service_sum;
		removed_weight_sum += removed_stats.removed_weight_sum;
		removed_volume_sum += removed_stats.removed_volume_sum;
		removed_depot_tt_sum += removed_stats.removed_depot_tt_sum;
		removed_position_sum += removed_stats.removed_position_sum;
	}

	for (const auto& r : stats.route_stats)
	{
		stats.total_score += r.collected_score;
		stats.total_served_customers += r.served_customers;

		sum_route_duration += r.route_duration;
		sum_route_time_util += r.route_time_utilization;
		sum_weight_util += r.weight_utilization;
		sum_volume_util += r.volume_utilization;
		sum_wait_per_customer += r.avg_waiting_time_per_customer;

		if (r.break_start_time >= 0.0)
		{
			sum_break_start += r.break_start_time;
			routes_with_break_start++;
		}
		if (r.break_position_norm >= 0.0)
		{
			sum_break_pos += r.break_position_norm;
			routes_with_break_pos++;
		}
		if (r.break_at_end_depot)
		{
			routes_break_at_end++;
		}
	}

	const int nRoutes = static_cast<int>(stats.route_stats.size());

	if (stats.total_served_customers > 0)
	{
		stats.avg_score_per_served_customer =
			static_cast<double>(stats.total_score) / stats.total_served_customers;
	}

	if (nRoutes > 0)
	{
		stats.avg_route_duration = sum_route_duration / nRoutes;
		stats.avg_route_time_utilization = sum_route_time_util / nRoutes;
		stats.avg_weight_utilization = sum_weight_util / nRoutes;
		stats.avg_volume_utilization = sum_volume_util / nRoutes;
		stats.avg_waiting_time_per_customer = sum_wait_per_customer / nRoutes;
		stats.pct_break_at_end_depot = 100.0 * static_cast<double>(routes_break_at_end) / nRoutes;
	}

	if (routes_with_break_start > 0)
		stats.avg_break_start_time = sum_break_start / routes_with_break_start;

	if (routes_with_break_pos > 0)
		stats.avg_break_position_norm = sum_break_pos / routes_with_break_pos;

	if (stats.total_removed_customers > 0)
	{
		stats.avg_removed_score =
			static_cast<double>(removed_score_sum) / stats.total_removed_customers;
		stats.avg_removed_tw_width =
			removed_tw_width_sum / stats.total_removed_customers;
		stats.avg_removed_service =
			removed_service_sum / stats.total_removed_customers;
		stats.avg_removed_weight =
			removed_weight_sum / stats.total_removed_customers;
		stats.avg_removed_volume =
			removed_volume_sum / stats.total_removed_customers;
		stats.avg_removed_depot_tt =
			removed_depot_tt_sum / stats.total_removed_customers;
		stats.avg_removed_position =
			removed_position_sum / stats.total_removed_customers;
	}

	return stats;
}

SolutionPerfStats Sol::collect_solution_stats() const
{
	SolutionPerfStats stats{};
	stats.route_stats.reserve(tours.size());
	stats.route_count = static_cast<int>(tours.size());

	double sum_route_duration = 0.0;
	double sum_route_time_util = 0.0;
	double sum_weight_util = 0.0;
	double sum_volume_util = 0.0;
	double sum_waiting = 0.0;
	double sum_customer_tw_slack = 0.0;
	double sum_break_start = 0.0;
	double sum_break_position = 0.0;
	int routes_with_break_start = 0;
	int routes_with_break_position = 0;
	int routes_break_at_end = 0;
	int total_customer_tw_slack_count = 0;
	double min_customer_tw_slack = DBL_MAX;
	double min_max_shift = DBL_MAX;

	for (const Tour& tour : tours)
	{
		RoutePerfStats route = tour.collect_route_stats();
		stats.route_stats.push_back(route);
		stats.customer_tw_slack_count += route.customer_tw_slack_count;
		stats.max_shift_count += route.max_shift_count;

		stats.total_score += route.collected_score;
		stats.total_served_customers += route.served_customers;
		if (route.served_customers > 0)
			++stats.nonempty_route_count;
		sum_route_duration += route.route_duration;
		sum_route_time_util += route.route_time_utilization;
		sum_weight_util += route.weight_utilization;
		sum_volume_util += route.volume_utilization;
		sum_waiting += route.total_waiting_time;
		if (route.customer_tw_slack_count > 0)
		{
			sum_customer_tw_slack +=
				route.avg_customer_tw_slack * route.customer_tw_slack_count;
			total_customer_tw_slack_count += route.customer_tw_slack_count;
		}

		stats.max_route_time_utilization =
			std::max(stats.max_route_time_utilization, route.route_time_utilization);
		stats.max_weight_utilization =
			std::max(stats.max_weight_utilization, route.weight_utilization);
		stats.max_volume_utilization =
			std::max(stats.max_volume_utilization, route.volume_utilization);
		if (route.route_time_utilization >= 0.90)
			++stats.routes_time_util_ge_90pct;
		if (route.weight_utilization >= 0.90)
			++stats.routes_weight_util_ge_90pct;
		if (route.volume_utilization >= 0.90)
			++stats.routes_volume_util_ge_90pct;
		if (route.customer_tw_slack_count > 0)
			min_customer_tw_slack =
			std::min(min_customer_tw_slack, route.min_customer_tw_slack);
		if (route.max_shift_count > 0)
			min_max_shift = std::min(min_max_shift, route.min_max_shift);

		if (route.break_start_time >= 0.0)
		{
			sum_break_start += route.break_start_time;
			++routes_with_break_start;
		}
		if (route.break_position_norm >= 0.0)
		{
			sum_break_position += route.break_position_norm;
			++routes_with_break_position;
		}
		if (route.break_at_end_depot)
			++routes_break_at_end;
	}

	const int route_count = static_cast<int>(stats.route_stats.size());
	if (stats.total_served_customers > 0)
	{
		stats.avg_score_per_served_customer =
			static_cast<double>(stats.total_score) / stats.total_served_customers;
		stats.avg_waiting_time_per_customer =
			sum_waiting / stats.total_served_customers;
	}
	stats.min_customer_tw_slack = min_customer_tw_slack < DBL_MAX
		? min_customer_tw_slack : -1.0;
	stats.min_max_shift = min_max_shift < DBL_MAX ? min_max_shift : -1.0;
	stats.avg_customer_tw_slack = total_customer_tw_slack_count > 0
		? sum_customer_tw_slack / total_customer_tw_slack_count : -1.0;

	if (route_count > 0)
	{
		// These means include empty route opportunities, so they describe average
		// utilization across the entire predefined fleet/day set.
		stats.avg_route_duration = sum_route_duration / route_count;
		stats.avg_route_time_utilization = sum_route_time_util / route_count;
		stats.avg_weight_utilization = sum_weight_util / route_count;
		stats.avg_volume_utilization = sum_volume_util / route_count;
		stats.pct_break_at_end_depot =
			100.0 * routes_break_at_end / route_count;
		stats.pct_routes_time_util_ge_90pct =
			100.0 * stats.routes_time_util_ge_90pct / route_count;
		stats.pct_routes_weight_util_ge_90pct =
			100.0 * stats.routes_weight_util_ge_90pct / route_count;
		stats.pct_routes_volume_util_ge_90pct =
			100.0 * stats.routes_volume_util_ge_90pct / route_count;
	}

	if (routes_with_break_start > 0)
		stats.avg_break_start_time = sum_break_start / routes_with_break_start;
	if (routes_with_break_position > 0)
		stats.avg_break_position_norm =
		sum_break_position / routes_with_break_position;

	return stats;
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

bool Tour::update_break()  // reposition break, update times and max_shift
{
	const int n = (int)seq.size();
	bool feasible = true;
	// reset actions
	std::fill(action.begin(), action.end(), 0);
	breakindex = -1;

	// --- forward pass WITHOUT break: arr0 (arrival-before-wait), svc (start-of-service) ---
	std::vector<double> arr0(n, 0.0);
	std::vector<double> svc(n, 0.0);

	double t0 = ins->t[index].EDT;
	deptime[0] = 0.0;

	for (int u = 0; u < n - 1; ++u) 
	{
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];

		double a = ins->arrival_time(o->con[p->index], t0); // arrival BEFORE waiting
		arr0[u + 1] = a;

		double s = (a < p->LTW[index] ? p->LTW[index] : a); // apply LTW wait
		svc[u + 1] = s;

		double leave = s + p->serv;                         // (no break in this phase)
		deptime[u + 1] = leave - ins->t[index].EDT;
		t0 = leave;
	}

	// compute max_shift on the no-break schedule
	calc_maxshift();

	// --- usable waiting AFTER breakstart at each vertex ---
	//     wait_after_start[j] = max(0, svc[j] - max(arr0[j], breakstart))
	std::vector<double> wait_after_start(n, 0.0);
	for (int j = 1; j < n; ++j) 
	{
		double usable = svc[j] - std::max(arr0[j], ins->breakstart);
		if (usable < 0.0) usable = 0.0;
		wait_after_start[j] = usable;
	}

	// --- earliest feasible index (start-only semantics) ---
	//     use arr0 (not current clock); ignore breakend here (we'll handle feasibility via slack)
	int earliestbreakindex = n - 1; // default to end depot
	for (int u = 0; u < n - 1; ++u) 
	{
		// Can begin by breakstart if we arrive after breakstart or can wait to LTW so that LTW-B >= breakstart
		double lhs = std::max(arr0[u + 1], seq[u + 1]->LTW[index] - ins->breakdur);
		if (lhs >= ins->breakstart - 1e-9) {
			earliestbreakindex = u + 1;
			break;
		}
	}

	// --- choose first position with enough budget; else force end depot (cannot fail) ---
	const double B = ins->breakdur;
	bool placed = false;
	const double eps = 1e-9;

	for (int u = earliestbreakindex; u < n; ++u)
	{
		const bool at_end_depot =(seq[u]->index == ins->maxvertices - 1);

		double break_start = arr0[u];
		if (!at_end_depot && break_start < ins->breakstart)
			break_start = ins->breakstart;

		if (break_start > ins->breakend + eps)
			continue;

		const double service_start =std::max(break_start + B, seq[u]->LTW[index]);

		if (service_start > seq[u]->UTW[index] + eps)
			continue;

		if (max_shift[u] + wait_after_start[u] + eps >= B)
		{
			action[u] = 1;
			breakindex = u;
			placed = true;
			break;
		}
	}
	if (!placed) 
	{
		//route is infeasible as break can not be placed
		feasible = false;
		return feasible;
		// unconditional fallback: end depot (special-case semantics will apply)
		//action.back() = 1;
		//breakindex = n - 1;
		//cout << "no feasible break position found, forcing end depot, solution might be infeasbible" << endl;
	}

	// --- forward recompute from breakindex-1 with break semantics ---
	double currenttime = (breakindex > 0)
		? deptime[breakindex - 1] + ins->t[index].EDT
		: ins->t[index].EDT;

	for (int u = std::max(0, breakindex - 1); u < n - 1; ++u) 
	{
		Ins::Vertex* o = seq[u];
		Ins::Vertex* p = seq[u + 1];

		double arr = ins->arrival_time(o->con[p->index], currenttime);

		if (action[u + 1]) 
		{
			const bool endp = (p->index == ins->maxvertices - 1);
			// start-only rule:
			//  - Non-depot: if early, wait to breakstart, then add B
			//  - Depot: allowed to start at arrival even if < breakstart
			if (!endp && arr < ins->breakstart) arr = ins->breakstart;
			arr += B; // take the break
		}

		// LTW after break if needed, then service
		if (arr < p->LTW[index]) arr = p->LTW[index];
		arr += p->serv;

		// local delta to max_shift (NO CLAMP; negative signals infeasibility)
		max_shift[u + 1] = (deptime[u + 1] + max_shift[u + 1]) - (arr - ins->t[index].EDT);

		deptime[u + 1] = arr - ins->t[index].EDT;
		currenttime = arr;
	}

	// --- safety: recompute end-depot from predecessor mirroring semantics ---
	{
		double tp  = ins->t[index].EDT + deptime[n - 2];
		double arr = ins->arrival_time(seq[n - 2]->con[seq[n - 1]->index], tp);
		if (action[n - 1]) arr += B;                             // depot break on arrival
		if (arr < seq[n - 1]->LTW[index]) arr = seq[n - 1]->LTW[index]; // honor LTW at depot (usually 0)
		arr += seq[n - 1]->serv;                                 // usually 0
		deptime[n - 1] = arr - ins->t[index].EDT;
	}

	// --- backward slack cap (unchanged in spirit) ---
	if (breakindex == n - 1) 
	{
		// break at end depot
		max_shift.back() = (ins->t[index].T_max - deptime.back());
		double arrive_cap = ins->t[index].LAT - B;
		if (ins->t[index].EDT + ins->t[index].T_max > ins->breakend + B) 
		{
			// If breakend is tighter, reflect that in the cap (still non-failing)
			max_shift.back() = (ins->breakend + B) - (deptime.back() + ins->t[index].EDT);
			arrive_cap = ins->breakend;
		}
		update_maxshift(0, n - 2, arrive_cap);
	} else 
	{
		// break before depot → cap by the successor’s latest arrival-before-service
		double cap_at_next = (ins->t[index].EDT + deptime[breakindex + 1] + max_shift[breakindex + 1]) - (seq[breakindex + 1]->serv + (action[breakindex + 1] ? B : 0.0));
		update_maxshift(0, breakindex, cap_at_next);
	}
	return feasible;
}

void Tour::update_maxshift(int start, int end, double arrivaltime)
{
	const double EDT = ins->t[index].EDT;
	const double B = ins->breakdur;
	const double BEND = ins->breakend;

	// Latest feasible absolute arrival at seq[end + 1],
	// before its break, waiting and service.
	double arr_cap = arrivaltime;

	for (int i = end; i > start; --i)
	{
		Ins::Vertex* y = seq[i];
		Ins::Vertex* z = seq[i + 1];
		const bool has_break_here = (action[i] == 1);

		// Latest departure from y that reaches z by arr_cap.
		double L_dep_i =
			ins->departure_time(y->con[z->index], arr_cap);

		// Service at y must start within its time window.
		L_dep_i = std::min(
			L_dep_i,
			y->UTW[index] + y->serv
		);

		// Departure slack, after break, waiting and service at y.
		max_shift[i] = L_dep_i - (EDT + deptime[i]);

		// Latest arrival at y, used for its predecessor.
		arr_cap = L_dep_i - y->serv;

		if (has_break_here)
		{
			arr_cap = std::min(arr_cap - B, BEND);
		}
	}
}

void Tour::calc_maxshift() 
{
	const int n = (int)seq.size();
	const double EDT = ins->t[index].EDT;          // == LTW at end depot
	const double LAT = ins->t[index].LAT;          // == EDT + T_max
	const double B = ins->breakdur;
	const double BEND = ins->breakend;

	// --- End depot slack on DEPARTURE time (finish time) ---
	// Latest allowed departure at end depot = min(LAT, BEND + (break?B:0))
	double L_dep_end = LAT;
	if (action.back() == 1) L_dep_end = std::min(L_dep_end, BEND + B);

	// deptime[] is relative to EDT, so absolute current departure is EDT + deptime.back()
	max_shift.back() = L_dep_end - (EDT + deptime.back());

	// --- Backward seed: latest ARRIVAL-BEFORE-SERVICE at end depot ---
	// If the break is at the depot, arrival must be ≤ min(LAT, BEND + B) - B
	double arr_cap = (action.back() == 1) ? (std::min(LAT, BEND + B) - B) : LAT;

	// --- Backward recursion over internal vertices ---
	for (int i = n - 2; i > 0; --i) 
	{
		Ins::Vertex* y = seq[i];
		Ins::Vertex* z = seq[i + 1];
		const bool has_break_here = (action[i] == 1);

		// Latest departure permitted by the successor and y's service window.
		double L_dep_i =ins->departure_time(y->con[z->index], arr_cap);

		L_dep_i = std::min(L_dep_i, y->UTW[index] + y->serv);

		// Departure slack: service and any break at y are already completed.
		max_shift[i] = L_dep_i - (EDT + deptime[i]);

		// Latest arrival at y, propagated backward to its predecessor.
		arr_cap = L_dep_i - y->serv;

		if (has_break_here)
		{
			arr_cap = std::min(arr_cap - B, BEND);
		}

	}
}

// Insert a previously validated customer after the given route position.
// Preconditions: candidate is unselected, position is before the end depot,
// and the original route's departures and slack cache are current.
void Tour::insert_vertex(Ins::Vertex* candidate, int position)
{
	const int inserted = position + 1;

	seq.insert(seq.begin() + inserted, candidate);
	deptime.insert(deptime.begin() + inserted, 0.0);
	action.insert(action.begin() + inserted, 0);
	max_shift.insert(max_shift.begin() + inserted, 0.0);

	score += candidate->score;
	volume += candidate->volume;
	weight += candidate->weight;
	if (position < breakindex)
		++breakindex;

	update(position, static_cast<int>(seq.size()));

	// The first OLD successor retains its suffix, so the delta-updated
	// departure slack at q gives a valid latest departure bound.
	const int q = inserted + 1;
	const double latest_dep_q =	ins->t[index].EDT + deptime[q] + max_shift[q];

	double latest_arrival_q = latest_dep_q - seq[q]->serv;
	if (action[q] == 1)
	{
		latest_arrival_q = std::min(
			latest_arrival_q - ins->breakdur,
			ins->breakend
		);
	}

	// Rebuild the inserted vertex AND all preceding customer positions.
	update_maxshift(0, inserted, latest_arrival_q);
}


bool Tour::remove_vertex(int position)
{
	bool feasible = true;
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
		feasible=update_break();
	}
	else
	{
		update(position - 1, int(seq.size()));//update travel time for all vertices after deletion
		calc_maxshift();
	}
	return feasible;
}//end remove_vertex tour version

bool Tour::remove_vertices(int position1, int position2)
{
	bool feasible = true;
	if (position1 == position2)
	{
		return feasible;
	}
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
	feasible = update_break();
	return feasible;
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

vector<double> Tour::compute_wait_suffix() const
{
	const int n = (int)seq.size();
	std::vector<double> arr0(n, 0.0), svc(n, 0.0), wait_suffix(n, 0.0);

	const double EDT = ins->t[index].EDT;
	const double BSTART = ins->breakstart;

	// --- no-break forward pass: arr0[j], svc[j]
	double t = deptime[0] + EDT;             // depart from start depot at scheduled time
	for (int j = 0; j < n - 1; ++j) 
	{
		auto* x = seq[j];
		auto* z = seq[j + 1];

		double a = ins->arrival_time(x->con[z->index], t); // arrival BEFORE any wait
		arr0[j + 1] = a;

		double s = (a < z->LTW[index] ? z->LTW[index] : a); // start-of-service with LTW (no break)
		svc[j + 1] = s;

		t = s + z->serv;                                    // depart after service (no break)
	}

	// --- usable waiting after breakstart, then suffix sum
	double acc = 0.0;
	for (int j = n - 1; j >= 1; --j) 
	{
		double t0 = std::max(arr0[j], BSTART);
		double usable = svc[j] - t0;
		if (usable < 0.0) usable = 0.0;
		acc += usable;
		wait_suffix[j] = acc;
	}
	return wait_suffix;
}

bool Tour::check()
{
	constexpr double EPS = 1e-6;
	bool tourok = true;

	const auto fail = [&](const char* message, int position = -1)
		{
			std::cout << term::fg(term::Color::red)<< "tour: " << index << " FAILURE: " << message;
			if (position >= 0)
				std::cout << " at position " << position;
			std::cout << '\n';
			tourok = false;
		};

	const auto mismatch = [&](const char* name, int position,
		double calculated, double stored)
		{
			std::cout << term::fg(term::Color::red)<< "tour: " << index << " FAILURE: " << name;
			if (position >= 0)
				std::cout << " at position " << position;
			std::cout << " calculated: " << calculated
				<< " stored: " << stored << '\n';
			tourok = false;
		};

	const int n = static_cast<int>(seq.size());
	if (n < 2)
	{
		fail("route must contain a start and an end depot");
		return false;
	}
	if (action.size() != seq.size() ||
		deptime.size() != seq.size() ||
		max_shift.size() != seq.size())
	{
		fail("seq/action/deptime/max_shift sizes do not match");
		return false;
	}
	for (int p = 0; p < n; ++p)
	{
		if (seq[p] == nullptr)
		{
			fail("null vertex pointer", p);
			return false;
		}
		if (!std::isfinite(deptime[p]))
		{
			fail("non-finite stored departure time", p);
			return false;
		}
		if (action[p] != 0 && action[p] != 1)
		{
			fail("break indicator must be zero or one", p);
			return false;
		}
	}

	const int end = n - 1;
	if (seq.front() != &ins->v[0] ||
		seq.back() != &ins->v[ins->maxvertices - 1])
	{
		fail("incorrect start or end depot");
		return false;
	}
	if (action[0] != 0)
	{
		fail("a break cannot be assigned to the start depot", 0);
		return false;
	}

	int amountbreaks = 0;
	int actual_break_position = -1;
	for (int p = 1; p <= end; ++p)
	{
		if (action[p] == 1)
		{
			++amountbreaks;
			actual_break_position = p;
		}
	}
	if (amountbreaks != 1)
	{
		fail("route must contain exactly one break");
		return false;
	}
	if (breakindex != actual_break_position)
		fail("breakindex does not match action[]");

	const auto& Td = ins->t[index];
	const double EDT = Td.EDT;
	const double B = ins->breakdur;
	const double BSTART = ins->breakstart;
	const double BEND = ins->breakend;

	// Recompute customer loads; check duplicate visits within this tour.
	double weightcheck = 0.0;
	double volumecheck = 0.0;
	std::vector<bool> seen(ins->maxvertices, false);
	for (int p = 1; p < end; ++p)
	{
		const int vertex = seq[p]->index;
		if (vertex <= 0 || vertex >= ins->maxvertices - 1)
		{
			fail("internal route position is not a customer", p);
			return false;
		}
		if (seen[vertex])
		{
			fail("customer occurs more than once", p);
			return false;
		}
		seen[vertex] = true;
		weightcheck += seq[p]->weight;
		volumecheck += seq[p]->volume;
	}
	if (!std::isfinite(weightcheck) || !std::isfinite(weight) ||
		std::abs(weightcheck - weight) > EPS)
		mismatch("weight", -1, weightcheck, weight);
	if (weightcheck > Td.W_max + EPS)
		fail("weight capacity exceeded");
	if (!std::isfinite(volumecheck) || !std::isfinite(volume) ||
		std::abs(volumecheck - volume) > EPS)
		mismatch("volume", -1, volumecheck, volume);
	if (volumecheck > Td.V_max + EPS)
		fail("volume capacity exceeded");

	// Forward validation, independently of cached max_shift values.
	std::vector<double> departurecheck(n, 0.0); // Absolute times.
	departurecheck[0] = EDT + deptime[0];
	if (deptime[0] < -EPS)
		fail("departure precedes the route-start reference", 0);

	double currenttime = departurecheck[0];
	for (int p = 1; p <= end; ++p)
	{
		Ins::Vertex* last = seq[p - 1];
		Ins::Vertex* current = seq[p];

		const double arrival =
			ins->arrival_time(last->con[current->index], currenttime);
		if (!std::isfinite(arrival))
		{
			fail("non-finite calculated arrival time", p);
			return false;
		}

		double ready = arrival;
		if (action[p] == 1)
		{
			// Customers may wait for BSTART. The end-depot exception
			// permits starting the break immediately upon arrival.
			const double break_start = (p == end)
				? arrival : std::max(arrival, BSTART);
			if (break_start > BEND + EPS)
				fail("break starts after breakend", p);
			ready = break_start + B;
		}

		const double service_start =
			std::max(ready, current->LTW[index]);
		if (service_start > current->UTW[index] + EPS)
			fail("service starts after the customer/depot closes", p);

		currenttime = service_start + current->serv;
		if (!std::isfinite(currenttime))
		{
			fail("non-finite calculated departure time", p);
			return false;
		}
		departurecheck[p] = currenttime;

		const double elapsed_departure = currenttime - EDT;
		if (std::abs(elapsed_departure - deptime[p]) > EPS)
			mismatch("departure time", p, elapsed_departure, deptime[p]);
	}

	if (currenttime - EDT > Td.T_max + EPS)
		fail("route-duration limit exceeded");
	if (currenttime > Td.LAT + EPS)
		fail("latest route-completion time exceeded");

	// Recompute departure slack without modifying the stored cache.
	// Use the independently reconstructed departures above.
	std::vector<double> max_shiftcheck(n, 0.0);
	double L_dep_end = Td.LAT;
	if (action[end] == 1)
		L_dep_end = std::min(L_dep_end, BEND + B);
	max_shiftcheck[end] = L_dep_end - departurecheck[end];

	double arr_cap = L_dep_end - ((action[end] == 1) ? B : 0.0);
	for (int p = end - 1; p > 0; --p)
	{
		Ins::Vertex* y = seq[p];
		Ins::Vertex* z = seq[p + 1];
		double L_dep =
			ins->departure_time(y->con[z->index], arr_cap);
		if (!std::isfinite(L_dep))
		{
			fail("non-finite latest departure time", p);
			return false;
		}
		L_dep = std::min(L_dep, y->UTW[index] + y->serv);

		// No customer-break cap on the post-service departure time.
		max_shiftcheck[p] = L_dep - departurecheck[p];

		// Apply the latest-break-start bound to incoming arrival time.
		arr_cap = L_dep - y->serv;
		if (action[p] == 1)
			arr_cap = std::min(arr_cap - B, BEND);
	}

	// Position 0 is excluded: calc_maxshift/update_maxshift do not set it.
	for (int p = 1; p <= end; ++p)
	{
		if (!std::isfinite(max_shift[p]) ||
			std::abs(max_shift[p] - max_shiftcheck[p]) > EPS)
			mismatch("max_shift", p, max_shiftcheck[p], max_shift[p]);
		if (max_shiftcheck[p] < -EPS)
			fail("negative recomputed slack", p);
	}

	return tourok;
}


std::pair<int, int> Tour::repair(RemovedVertexStats& removed_stats)
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

		Ins::Vertex* removed_vertex = seq[end - 1];

		// log removed-customer characteristics
		removed_stats.removed_count++;
		removed_stats.removed_score_sum += removed_vertex->score;
		removed_stats.removed_service_sum += removed_vertex->serv;
		removed_stats.removed_weight_sum += removed_vertex->weight;
		removed_stats.removed_volume_sum += removed_vertex->volume;
		removed_stats.removed_tw_width_sum += (removed_vertex->UTW[index] - removed_vertex->LTW[index]);

		// time-independent depot travel time from start depot to removed customer
		Ins::Vertex* start_depot = seq.front();
		removed_stats.removed_depot_tt_sum += start_depot->con[removed_vertex->index]->determin;

		// normalized position in route among current regular vertices
		int served_before = std::max(1, static_cast<int>(seq.size()) - 2);
		double pos_norm = static_cast<double>(end - 1) / served_before;
		removed_stats.removed_position_sum += pos_norm;

		scoredecrease += removed_vertex->score;
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

RoutePerfStats Tour::collect_route_stats() const
{
	RoutePerfStats s{};
	// Set explicit sentinels here as well as in the struct definition. This
	// keeps the audit safe if these stats are declared without default values.
	s.min_customer_tw_slack = -1.0;
	s.avg_customer_tw_slack = -1.0;
	s.customer_tw_slack_count = 0;
	s.customers_tw_slack_le_0_1 = 0;
	s.min_max_shift = -1.0;
	s.max_shift_count = 0;

	s.route_index = index;
	s.collected_score = score;
	s.served_customers = std::max(0, static_cast<int>(seq.size()) - 2);
	const int n = static_cast<int>(seq.size());
	if (n < 2 || deptime.size() != seq.size() || action.size() != seq.size())
		throw std::logic_error("Cannot collect route audit: route schedule arrays have inconsistent sizes");

	s.route_time_budget = ins->t[index].T_max;

	const double EDT = ins->t[index].EDT;
	double currenttime = EDT + deptime[0];
	double total_waiting = 0.0;
	double break_start_abs = -1.0;
	double customer_tw_slack_sum = 0.0;
	double customer_tw_slack_min = DBL_MAX;
	int customer_tw_slack_count = 0;
	int customers_tw_slack_le_01 = 0;

	int end = n - 1;

	for (int i = 0; i < end; ++i)
	{
		Ins::Vertex* last = seq[i];
		Ins::Vertex* current = seq[i + 1];
		int breakcurrent = action[i + 1];

		const double raw_arrival =
			ins->arrival_time(last->con[current->index], currenttime);
		double service_start = raw_arrival;

		if (breakcurrent == 1)
		{
			const bool at_end_depot =
				current->index == ins->maxvertices - 1;
			double break_start = raw_arrival;
			if (!at_end_depot && break_start < ins->breakstart)
				break_start = ins->breakstart;

			if (break_start_abs < 0.0)
				break_start_abs = break_start;

			total_waiting += std::max(0.0, break_start - raw_arrival);
			const double break_end = break_start + ins->breakdur;
			service_start = std::max(break_end, current->LTW[index]);
			total_waiting += std::max(0.0, service_start - break_end);
		}
		else
		{
			service_start = std::max(raw_arrival, current->LTW[index]);
			total_waiting += std::max(0.0, service_start - raw_arrival);
		}

		// Record service-window slack for customer visits only (exclude both depots).
		if (i + 1 > 0 && i + 1 < end)
		{
			const double tw_slack = current->UTW[index] - service_start;
			customer_tw_slack_sum += tw_slack;
			customer_tw_slack_min = std::min(customer_tw_slack_min, tw_slack);
			++customer_tw_slack_count;
			if (tw_slack <= 0.1 + 1e-9)
				++customers_tw_slack_le_01;
		}

		currenttime = service_start + current->serv;
	}

	s.route_duration = currenttime - EDT;
	if (s.route_time_budget > 1e-9)
		s.route_time_utilization = s.route_duration / s.route_time_budget;

	s.weight_used = weight;
	s.weight_capacity = ins->t[index].W_max;
	s.weight_slack = s.weight_capacity - s.weight_used;
	if (s.weight_capacity > 1e-9)
		s.weight_utilization = s.weight_used / s.weight_capacity;

	s.volume_used = volume;
	s.volume_capacity = ins->t[index].V_max;
	s.volume_slack = s.volume_capacity - s.volume_used;
	if (s.volume_capacity > 1e-9)
		s.volume_utilization = s.volume_used / s.volume_capacity;

	s.min_customer_tw_slack = customer_tw_slack_count > 0
		? customer_tw_slack_min : -1.0;
	s.avg_customer_tw_slack = customer_tw_slack_count > 0
		? customer_tw_slack_sum / customer_tw_slack_count : -1.0;
	s.customer_tw_slack_count = customer_tw_slack_count;
	s.customers_tw_slack_le_0_1 = customers_tw_slack_le_01;

	// Recompute a local copy of the backward slack values. This avoids reporting
	// stale or never-updated max_shift entries while leaving the solution intact.
	if (max_shift.size() != seq.size())
		throw std::logic_error("Cannot collect route audit: max_shift size does not match the route");
	Tour slack_copy = *this;
	slack_copy.calc_maxshift();
	double min_shift = DBL_MAX;
	int max_shift_count = 0;
	for (int p = 1; p < n; ++p) // all visited vertices after the start depot, including end depot
	{
		const double slack = slack_copy.max_shift[p];
		if (std::isfinite(slack))
		{
			min_shift = std::min(min_shift, slack);
			++max_shift_count;
		}
	}
	s.min_max_shift = min_shift < DBL_MAX ? min_shift : -1.0;
	s.max_shift_count = max_shift_count;

	s.total_waiting_time = total_waiting;
	if (s.served_customers > 0)
		s.avg_waiting_time_per_customer = total_waiting / s.served_customers;

	s.break_taken = (breakindex >= 0 && breakindex < static_cast<int>(seq.size()));
	s.break_at_end_depot = (breakindex == static_cast<int>(seq.size()) - 1);

	if (break_start_abs >= 0.0)
	{
		s.break_start_time = break_start_abs - ins->t[index].EDT;
		s.break_start_slack = ins->breakend - break_start_abs;
	}

	if (s.served_customers > 0 && s.break_taken)
	{
		// normalize by number of customer visits
		if (breakindex <= 0)
			s.break_position_norm = 0.0;
		else if (breakindex >= static_cast<int>(seq.size()) - 1)
			s.break_position_norm = 1.0;
		else
			s.break_position_norm = static_cast<double>(breakindex) / s.served_customers;
	}

	return s;
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

bool Sol::remove_vertex(Sol::Tour& tour, int position)
{
	bool feasible = true;
	Ins::Vertex* candidate = tour.seq[position];
	available[candidate->index] = true;//make vertex back available
	score -= candidate->score;// update score of the new solution
	feasible = tour.remove_vertex(position);
	return feasible;
}

bool Sol::remove_vertices(Tour& tour, int position1, int position2)
{
	bool feasible = true;
	Ins::Vertex* candidate1 = tour.seq[position1];
	Ins::Vertex* candidate2 = tour.seq[position2];
	available[candidate1->index] = true;//make vertex back available
	score -= candidate1->score;// update score of the new solution
	available[candidate2->index] = true;//make vertex back available
	score -= candidate2->score;// update score of the new solution
	feasible = tour.remove_vertices(position1,position2);
	return feasible;
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
