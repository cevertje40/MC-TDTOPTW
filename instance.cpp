#include "instance.h"

using namespace std;

void Ins::read_time_independent_traveltime()
{
	// read in time-independent (freeflow) travel time
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	ifstream titt;
	titt.open(filepath + "titt" + to_string(maxvertices) + ".TXT", ifstream::in);
	if (titt.is_open())
	{
		//cout << "reading time-independent travel time" << endl;
		//assign connections to vertex objects and read free flow travel time
		for (int i = 0; i < maxvertices; ++i)
		{
			v[i].con.resize(maxvertices);
			string line;
			getline(titt, line);
			stringstream str(line);
			for (int j = 0; j < maxvertices; ++j)
			{
				int counter = i * maxvertices + j;
				c[counter].from = i;
				c[counter].to = j;
				str >> c[counter].determin;
				str.ignore();
				v[i].con[j] = &c[counter];
			}
		}
		titt.close();
	}
	else
	{
		cout << term::fg(term::Color::red) << "could not open time-independent travel time file" << endl;
	}
}//end read time independent


void Ins::read_time_dependent_traveltime()
{
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	ifstream tt(filepath + "tt" + to_string(maxvertices) + ".TXT");
	if (!tt.is_open()) {
		cout << term::fg(term::Color::red) << "could not open time-dependent travel time file" << endl;
		return;
	}
	//cout << "reading time-dependent travel time" << endl;
	for (int i = 0; i < maxvertices; ++i) {
		for (int j = 0; j < maxvertices; ++j) {
			string line;
			getline(tt, line);
			stringstream str(line);

			vector<double> dump(maxtimeslots);
			for (double& val : dump) {
				str >> val;
				str.ignore();  // ignore delimiter (assumes comma or space)
			}

			auto& mu = v[i].con[j]->mu;
			auto& oneplusmu = v[i].con[j]->oneplusmu;
			auto& nu = v[i].con[j]->nu;
			oneplusmu.resize(maxtimeslots);
			mu.resize(maxtimeslots);
			nu.resize(maxtimeslots);
			//double mintraveltime = DBL_MAX;
			for (int t = 0; t < maxtimeslots; ++t) 
			{
				//mintraveltime = min(mintraveltime, dump[t]);
				if (t == maxtimeslots - 1) {
					mu[t] = 0.0;
					oneplusmu[t] = 1.0;
					nu[t] = dump[t];
				}
				else {
					double delta_time = time_periods[t + 1] - time_periods[t];
					mu[t] = (dump[t + 1] - dump[t]) / delta_time;
					oneplusmu[t] = 1.0 + mu[t];
					nu[t] = dump[t] - mu[t] * time_periods[t];
				}
			}
			//v[i].con[j]->determin = mintraveltime;//temporary test
		}
	}
	tt.close();
}//end read time dependent travel time


Ins::Ins(MCTDTOPTW textfile)
{
	//read in vertex and tour information from txt file and populate v and t objects
	ifstream ifs;
	ifs.open(textfile.path + textfile.name, ifstream::in);
	if (ifs.is_open())
	{
		string line;
		getline(ifs, line);
		stringstream str(line);//first line
		str >> maxvertices;//maxvertices
		getline(ifs, line);
		str = stringstream(line);//second line
		str >> maxtours;//maxtours
		v.resize(maxvertices);
		t.resize(maxtours);
		c.resize(maxvertices * maxvertices);
		getline(ifs, line);//third line, break info
		str = stringstream(line);
		str >> breakdur >> breakstart >> breakend;
		maxscore = 0;
		for (int vertex = 0; vertex < maxvertices; ++vertex)
		{
			getline(ifs, line);//maxvertices lines vertex info
			str = stringstream(line);
			v[vertex].index = vertex;
			str >> v[vertex].id >> v[vertex].score >> v[vertex].serv >> v[vertex].weight >> v[vertex].volume;
			v[vertex].serv = std::round(v[vertex].serv * 100.0) / 100.0;
			v[vertex].weight = std::round(v[vertex].weight * 100.0) / 100.0;
			v[vertex].volume = std::round(v[vertex].volume * 100.0) / 100.0;
			maxscore += v[vertex].score;
			v[vertex].LTW.resize(maxtours);
			v[vertex].UTW.resize(maxtours);
			for (int tour = 0; tour < maxtours; ++tour)
			{
				str >> v[vertex].LTW[tour] >> v[vertex].UTW[tour];
			}
		}
		for (int tour = 0; tour < maxtours; ++tour)
		{
			getline(ifs, line);//maxtours lines tour info
			str = stringstream(line);
			str >> t[tour].T_max >> t[tour].W_max >> t[tour].V_max;
			t[tour].index = tour;
			t[tour].id = tour + 1;
			t[tour].EDT = v.back().LTW[tour];
			t[tour].LAT = t[tour].EDT + t[tour].T_max;
			t[tour].startv = &v[0];
			t[tour].endv = &v[maxvertices - 1];
		}
		ifs.close();
	}
	else
	{
		cout << " can not open MC-TDTOPTW instance file" << endl;
	}
}

Ins::Ins(CTOP textfile)
{
	ifstream ifs;
	string filepath = textfile.path + textfile.name;
	ifs.open(filepath, ifstream::in);
	if (ifs.is_open())
	{
		string line;
		string dump;
		vector<double> x_coordinates;//initialised in read_in
		vector<double> y_coordinates;//initialised in read_in
		double C_max;
		double T_max;
		getline(ifs, line);//first line contains name
		getline(ifs, line);//second line is empty
		getline(ifs, line);//third line max vehicles
		stringstream str(line);//store line as stringstream
		str >> dump;
		str >> maxtours;
		getline(ifs, line);//fourth line maximum capacity
		str = stringstream(line);
		str >> dump;
		str >> C_max;
		getline(ifs, line);//fifth line maximum time
		str = stringstream(line);//store line as stringstream
		str >> dump;
		str >> T_max;
		getline(ifs, line);//sixth line is empty
		getline(ifs, line);//seventh line depot coordinates
		str = stringstream(line);//store line as stringstream
		str >> dump;
		double depotx = 0.0;
		double depoty = 0.0;
		str >> depotx;
		str >> depoty;
		getline(ifs, line);//eight line is empty
		getline(ifs, line);//ninth line is maxvertices
		str = stringstream(line);//store line as stringstream
		str >> dump;
		str >> maxvertices;
		maxvertices += 2;//add two vertices for the depots
		x_coordinates.resize(maxvertices, 0);
		y_coordinates.resize(maxvertices, 0);
		x_coordinates[0] = depotx;
		y_coordinates[0] = depoty;
		x_coordinates[maxvertices - 1] = depotx;
		y_coordinates[maxvertices - 1] = depoty;
		v.resize(maxvertices);
		t.resize(maxtours);
		c.resize(maxvertices * maxvertices);
		getline(ifs, line);//tenth line is empty
		getline(ifs, line);//eleventh line is text
		//regular vertices creation
		maxscore = 0;
		for (int i = 1; i < maxvertices - 1; ++i)//read maxvertices-2 amount regular vertices
		{
			v[i].id = i;
			v[i].index = i;
			getline(ifs, line);
			str = stringstream(line);
			str >> x_coordinates[i];
			str >> y_coordinates[i];
			str >> v[i].weight;
			v[i].volume = v[i].weight;
			str >> v[i].serv;
			str >> v[i].score;
			maxscore += v[i].score;
			//v[i].serv = 0.0;//set service time equal to 0 for set1-3, comment out for set4-6
			v[i].LTW.resize(maxtours);
			v[i].UTW.resize(maxtours);
			for (int tour = 0; tour < maxtours; ++tour)
			{
				v[i].LTW[tour] = time_periods[0];
				v[i].UTW[tour] = time_periods[0] + T_max;
			}
		}
		ifs.close();
		//start & end depot creation
		v[0].id = 0;
		v[0].index = 0;
		v[0].score = 0;
		v[0].serv = 0;
		v[0].weight = 0;
		v[0].volume = 0;
		v[0].LTW.resize(maxtours);
		v[0].UTW.resize(maxtours);
		for (int tour = 0; tour < maxtours; ++tour)
		{
			v[0].LTW[tour] = time_periods[0];
			v[0].UTW[tour] = time_periods[0] + T_max;
		}
		v[maxvertices - 1].id = maxvertices - 1;
		v[maxvertices - 1].index = maxvertices - 1;
		v[maxvertices - 1].score = 0;
		v[maxvertices - 1].serv = 0;
		v[maxvertices - 1].weight = 0;
		v[maxvertices - 1].volume = 0;
		v[maxvertices - 1].LTW.resize(maxtours);
		v[maxvertices - 1].UTW.resize(maxtours);
		for (int tour = 0; tour < maxtours; ++tour)
		{
			v[maxvertices - 1].LTW[tour] = time_periods[0];
			v[maxvertices - 1].UTW[tour] = time_periods[0] + T_max;
		}

		// route setup
		for (int tour = 0; tour < maxtours; ++tour)
		{
			t[tour].index = tour;
			t[tour].id = tour;
			t[tour].T_max = T_max;
			t[tour].W_max = C_max;//CTOP only has one capacity constraint
			t[tour].V_max = C_max;//CTOP only has one capacity constraint
			t[tour].EDT = time_periods[0];
			t[tour].LAT = time_periods[0] + T_max;
			t[tour].startv = &v[0];
			t[tour].endv = &v[maxvertices - 1];
		}
		//break setup
		breakstart = time_periods[0];
		breakdur = 0.0;
		breakend = DBL_MAX;

		//construct travel time matrix
		for (int i = 0; i < maxvertices; ++i)
		{
			v[i].con.resize(maxvertices);
			for (int j = 0; j < maxvertices; ++j)
			{
				int counter = i * maxvertices + j;
				c[counter].from = i;
				c[counter].to = j;
				//calculate EUC_2D: Euclidean distance rounded
				double dx = x_coordinates[i] - x_coordinates[j];
				double dy = y_coordinates[i] - y_coordinates[j];
				c[counter].determin = static_cast<int>(std::floor(std::hypot(dx, dy) + 0.5));
				c[counter].mu.resize(maxtimeslots);
				c[counter].nu.resize(maxtimeslots);
				c[counter].oneplusmu.resize(maxtimeslots);
				for (int t = 0; t < maxtimeslots; ++t)
				{
					c[counter].mu[t] = 0.0;
					c[counter].oneplusmu[t] = 1.0;
					c[counter].nu[t] = c[counter].determin;
				}
				v[i].con[j] = &c[counter];
			}
		}
	}
	else
	{
		cout << endl << " can not open CTOP instance file" << endl;;
	}
}

struct PlantedReplayResult
{
	double completion_time = 0.0;
	double loadW = 0.0;
	double loadV = 0.0;
	bool break_taken = false;
	bool feasible = true;
};

Ins::PlantedReplayResult Ins::replay_planted_route_break_first_customer(const std::vector<int>& route,int tour_idx,bool verbose) const
{
	PlantedReplayResult res;
	double currenttime = t[tour_idx].EDT;

	if (verbose)
	{
		std::cout << "Checking planted route " << tour_idx << std::endl;
	}

	for (size_t k = 0; k + 1 < route.size(); ++k)
	{
		int source = route[k];
		int target = route[k + 1];

		int slot = find_t(currenttime);
		double arrivaltime = arrival_time(v[source].con[target], currenttime);

		if (verbose)
		{
			std::cout << "  arc " << source << " -> " << target
				<< " dep=" << currenttime
				<< " slot=" << slot
				<< " tt=" << (arrivaltime - currenttime)
				<< " arr=" << arrivaltime;
		}

		if (target != maxvertices - 1)
		{
			// fixed rule: break at the first customer only
			if (!res.break_taken && k == 0)
			{
				if (verbose)
					std::cout << " break@" << target << " start=" << arrivaltime;

				// optional strict check: break must lie in break window
				if (arrivaltime < breakstart || arrivaltime > breakend)
				{
					res.feasible = false;
					if (verbose)
						std::cout << " BREAK-WINDOW-FAIL";
				}

				arrivaltime += breakdur;
				res.break_taken = true;
			}

			if (arrivaltime < v[target].LTW[tour_idx])
			{
				if (verbose)
					std::cout << " wait=" << (v[target].LTW[tour_idx] - arrivaltime);
				arrivaltime = v[target].LTW[tour_idx];
			}

			if (arrivaltime > v[target].UTW[tour_idx] + 1e-9)
			{
				res.feasible = false;
				if (verbose)
					std::cout << " TW-FAIL";
			}

			res.loadW += v[target].weight;
			res.loadV += v[target].volume;

			arrivaltime += v[target].serv;

			if (verbose)
			{
				std::cout << " serv=" << v[target].serv
					<< " dep_next=" << arrivaltime;
			}
		}

		if (verbose)
			std::cout << std::endl;

		currenttime = arrivaltime;
	}

	res.completion_time = currenttime;

	// if route had at least one customer, break should have been taken
	if (route.size() > 2 && !res.break_taken)
		res.feasible = false;

	return res;
}

static std::vector<std::vector<char>> build_planted_arc_matrix(int maxvertices,const std::vector<std::vector<int>>& full_routes)
{
	std::vector<std::vector<char>> planted_arc(	maxvertices,std::vector<char>(maxvertices, 0));

	for (const auto& route : full_routes)
	{
		for (size_t k = 0; k + 1 < route.size(); ++k)
		{
			int i = route[k];
			int j = route[k + 1];
			planted_arc[i][j] = 1;
		}
	}

	return planted_arc;
}

Ins::Ins(KnownOptimalCTOP data)
{
	ifstream ifs;
	string filepath = data.path + data.name;
	ifs.open(filepath, ifstream::in);

	if (!ifs.is_open())
	{
		cout << endl << " can not open CTOP instance file" << endl;
		return;
	}

	string line;
	string dump;
	vector<double> x_coordinates;
	vector<double> y_coordinates;
	double C_max = 0.0;
	double T_max_ctop = 0.0;

	// ------------------------------------------------------------
	// 1. Read CTOP file
	// ------------------------------------------------------------
	getline(ifs, line);
	getline(ifs, line);

	getline(ifs, line);
	stringstream str(line);
	str >> dump;
	str >> maxtours;

	getline(ifs, line);
	str = stringstream(line);
	str >> dump;
	str >> C_max;

	getline(ifs, line);
	str = stringstream(line);
	str >> dump;
	str >> T_max_ctop;

	getline(ifs, line);
	getline(ifs, line);
	str = stringstream(line);
	str >> dump;
	double depotx = 0.0;
	double depoty = 0.0;
	str >> depotx;
	str >> depoty;

	getline(ifs, line);
	getline(ifs, line);
	str = stringstream(line);
	str >> dump;
	str >> maxvertices;
	maxvertices += 2;

	x_coordinates.resize(maxvertices, 0.0);
	y_coordinates.resize(maxvertices, 0.0);
	x_coordinates[0] = depotx;
	y_coordinates[0] = depoty;
	x_coordinates[maxvertices - 1] = depotx;
	y_coordinates[maxvertices - 1] = depoty;

	v.resize(maxvertices);
	t.resize(maxtours);
	c.resize(maxvertices * maxvertices);

	getline(ifs, line);
	getline(ifs, line);

	maxscore = 0.0;

	// ------------------------------------------------------------
	// 2. Create regular vertices
	// ------------------------------------------------------------
	for (int i = 1; i < maxvertices - 1; ++i)
	{
		v[i].id = i;
		v[i].index = i;

		getline(ifs, line);
		str = stringstream(line);

		str >> x_coordinates[i];
		str >> y_coordinates[i];
		str >> v[i].weight;
		v[i].volume = v[i].weight;
		str >> v[i].serv;
		str >> v[i].score;

		v[i].serv *= data.time_scale;

		maxscore += v[i].score;

		v[i].LTW.resize(maxtours, 0.0);
		v[i].UTW.resize(maxtours, 0.0);
	}
	ifs.close();

	// ------------------------------------------------------------
	// 3. Create depots
	// ------------------------------------------------------------
	v[0].id = 0;
	v[0].index = 0;
	v[0].score = 0;
	v[0].serv = 0.0;
	v[0].weight = 0.0;
	v[0].volume = 0.0;
	v[0].LTW.resize(maxtours, 0.0);
	v[0].UTW.resize(maxtours, 0.0);

	v[maxvertices - 1].id = maxvertices - 1;
	v[maxvertices - 1].index = maxvertices - 1;
	v[maxvertices - 1].score = 0;
	v[maxvertices - 1].serv = 0.0;
	v[maxvertices - 1].weight = 0.0;
	v[maxvertices - 1].volume = 0.0;
	v[maxvertices - 1].LTW.resize(maxtours, 0.0);
	v[maxvertices - 1].UTW.resize(maxtours, 0.0);

	// ------------------------------------------------------------
	// 4. Build deterministic scaled CTOP matrix
	// ------------------------------------------------------------
	for (int i = 0; i < maxvertices; ++i)
	{
		v[i].con.resize(maxvertices);

		for (int j = 0; j < maxvertices; ++j)
		{
			int counter = i * maxvertices + j;
			c[counter].from = i;
			c[counter].to = j;

			double dx = x_coordinates[i] - x_coordinates[j];
			double dy = y_coordinates[i] - y_coordinates[j];
			double euc = static_cast<int>(std::floor(std::hypot(dx, dy) + 0.5));

			c[counter].determin = euc * data.time_scale;

			c[counter].mu.resize(maxtimeslots, 0.0);
			c[counter].oneplusmu.resize(maxtimeslots, 1.0);
			c[counter].nu.resize(maxtimeslots, c[counter].determin);

			v[i].con[j] = &c[counter];
		}
	}

	// ------------------------------------------------------------
	// 5. Build full planted routes [0 ... end depot]
	// ------------------------------------------------------------
	vector<vector<int>> full_routes(maxtours);
	vector<vector<int>> inputb(maxtours);
	std::string name = "sol_" + data.name;
	
	ifs.open(name, ifstream::in);
	if (ifs.is_open())
	{
		for (int t = 0; t < maxtours; ++t)
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
				full_routes[t].push_back(index);
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
	
	// ------------------------------------------------------------
	// 6. Record route loads
	// ------------------------------------------------------------
	std::vector<double> planted_route_weight(full_routes.size(), 0.0);
	std::vector<double> planted_route_volume(full_routes.size(), 0.0);

	for (size_t r = 0; r < full_routes.size(); ++r)
	{
		for (size_t k = 1; k + 1 < full_routes[r].size(); ++k)
		{
			int cust = full_routes[r][k];
			planted_route_weight[r] += v[cust].weight;
			planted_route_volume[r] += v[cust].volume;
		}
	}

	// ------------------------------------------------------------
	// 7. Break settings
	// ------------------------------------------------------------
	breakdur = data.break_dur;
	breakstart = time_periods[0] + data.break_start;
	breakend = time_periods[0] + data.break_end;

	// ------------------------------------------------------------
	// 8. Common CTOP-like Tmax for all tours
	// ------------------------------------------------------------
	double common_Tmax = T_max_ctop * data.time_scale + breakdur + data.horizon_slack;

	// ------------------------------------------------------------
	// 9. Finalize tours immediately with common Tmax
	// ------------------------------------------------------------
	for (int tour = 0; tour < maxtours; ++tour)
	{
		t[tour].index = tour;
		t[tour].id = tour;
		t[tour].startv = &v[0];
		t[tour].endv = &v[maxvertices - 1];
		t[tour].EDT = time_periods[0];
		t[tour].T_max = common_Tmax;
		t[tour].LAT = t[tour].EDT + t[tour].T_max;

		if (tour < static_cast<int>(full_routes.size()))
		{
			t[tour].W_max = C_max;
			t[tour].V_max = C_max;
		}
		else
		{
			if (data.disable_unused_tours)
			{
				t[tour].W_max = 0.0;
				t[tour].V_max = 0.0;
			}
			else
			{
				t[tour].W_max = C_max;
				t[tour].V_max = C_max;
			}
		}
	}

	// ------------------------------------------------------------
	// 10. Wide-open time windows
	// ------------------------------------------------------------
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int tour = 0; tour < maxtours; ++tour)
		{
			v[i].LTW[tour] = t[tour].EDT;
			v[i].UTW[tour] = t[tour].LAT;
		}
	}

	// ------------------------------------------------------------
	// 11. Mark planted arcs
	// ------------------------------------------------------------
	std::vector<std::vector<char>> planted_arc =
		build_planted_arc_matrix(maxvertices, full_routes);

	// ------------------------------------------------------------
	// 12. Build globally affine arc profiles
	//     Same coefficients over all slots for a given arc
	// ------------------------------------------------------------
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int j = 0; j < maxvertices; ++j)
		{
			int counter = i * maxvertices + j;
			double base = c[counter].determin;

			double mu_arc = 0.0;
			double nu_arc = base;

			if (planted_arc[i][j])
			{
				mu_arc = data.planted_mu;
				nu_arc = base + data.planted_extra_nu;
			}
			else
			{
				mu_arc = data.nonplanted_mu;
				nu_arc = base + data.nonplanted_extra_nu_factor * base;
			}

			// safety: keep tt >= determin at route start
			double min_tau = time_periods[0];
			if (mu_arc * min_tau + nu_arc < base)
			{
				nu_arc = base - mu_arc * min_tau;
			}

			for (int ts = 0; ts < maxtimeslots; ++ts)
			{
				c[counter].mu[ts] = mu_arc;
				c[counter].oneplusmu[ts] = 1.0 + mu_arc;
				c[counter].nu[ts] = nu_arc;
			}
		}
	}

	// ------------------------------------------------------------
	// 13. Diagnostics: FIFO / continuity / lower bound
	// ------------------------------------------------------------
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int j = 0; j < maxvertices; ++j)
		{
			Ins::Connec* arc = v[i].con[j];
			double base = arc->determin;

			for (int ts = 0; ts < maxtimeslots; ++ts)
			{
				if (arc->oneplusmu[ts] < -1e-9)
				{
					std::cout << "FIFO slope violation on arc "
						<< i << " -> " << j
						<< " slot " << ts
						<< " oneplusmu = " << arc->oneplusmu[ts] << std::endl;
				}

				double tau0 = time_periods[ts];
				double tau1 = time_periods[ts + 1];

				double tt0 = arc->mu[ts] * tau0 + arc->nu[ts];
				double tt1 = arc->mu[ts] * tau1 + arc->nu[ts];

				if (tt0 < base - 1e-9)
				{
					std::cout << "Travel time below deterministic on arc "
						<< i << " -> " << j
						<< " slot " << ts
						<< " at left boundary: tt=" << tt0
						<< " determin=" << base << std::endl;
				}

				if (tt1 < base - 1e-9)
				{
					std::cout << "Travel time below deterministic on arc "
						<< i << " -> " << j
						<< " slot " << ts
						<< " at right boundary: tt=" << tt1
						<< " determin=" << base << std::endl;
				}
			}

			for (int ts = 0; ts < maxtimeslots - 1; ++ts)
			{
				double boundary = time_periods[ts + 1];
				double left = arc->oneplusmu[ts] * boundary + arc->nu[ts];
				double right = arc->oneplusmu[ts + 1] * boundary + arc->nu[ts + 1];

				if (right < left - 1e-9)
				{
					std::cout << "FIFO boundary violation on arc "
						<< i << " -> " << j
						<< " between slots " << ts
						<< " and " << ts + 1
						<< " left=" << left
						<< " right=" << right << std::endl;
				}
			}
		}
	}

	// ------------------------------------------------------------
	// 14. Planted-route diagnostics only
	//     Use your fixed first-customer break replay
	// ------------------------------------------------------------
	for (size_t r = 0; r < full_routes.size(); ++r)
	{
		PlantedReplayResult rr =
			replay_planted_route_break_first_customer(full_routes[r], (int)r, false);

		if (!rr.feasible)
		{
			std::cout << "Planted route infeasible under fixed first-customer break rule on route "
				<< r << std::endl;
		}

		if (rr.loadW > t[r].W_max + 1e-9)
		{
			std::cout << "Weight violation on planted route " << r << std::endl;
		}

		if (rr.loadV > t[r].V_max + 1e-9)
		{
			std::cout << "Volume violation on planted route " << r << std::endl;
		}

		double total = rr.completion_time - t[r].EDT;
		if (total > t[r].T_max + 1e-9)
		{
			std::cout << "Duration violation on planted route " << r
				<< " total=" << total
				<< " Tmax=" << t[r].T_max << std::endl;
		}
	}

	// ------------------------------------------------------------
	// 15. Basic diagnostics
	// ------------------------------------------------------------
	if (full_routes.size() > static_cast<size_t>(maxtours))
	{
		cout << "Warning: more planted routes than available tours." << endl;
	}

	for (size_t r = 0; r < full_routes.size(); ++r)
	{
		if (planted_route_weight[r] > C_max)
		{
			cout << "Warning: planted route " << r << " exceeds weight capacity." << endl;
		}
		if (planted_route_volume[r] > C_max)
		{
			cout << "Warning: planted route " << r << " exceeds volume capacity." << endl;
		}
	}
}

void Ins::construct_time_independent_traveltime(Graph& graph)
{
	clock_t start, end;
	start = clock();
	vector<int>mapper(maxvertices, 0);
	vector<int>destinations;
	vector<vector<double>> dump;
	for (int i = 0; i < maxvertices; ++i)
	{
		if (find(destinations.begin(),destinations.end(), v[i].id - 1) == destinations.end())//unique beindex
		{
			destinations.push_back(v[i].id - 1);//only add unique beindices to target list
			mapper[i] = int(destinations.size()) - 1;
		}
		else//non-unique beindex
		{
			//map this vertex to the first occurence in the desination vector
			mapper[i] = int(find(destinations.begin(), destinations.end(), v[i].id - 1)-destinations.begin());
		}
		dump.push_back(vector<double>());//1 target array per vertex
	}
	cout << "construct time independent travel time" << endl;
	#pragma omp parallel num_threads(12)
	{//start parallel session
		#pragma omp for nowait
		for (int i = 0; i < maxvertices; ++i)
		{
			int thread = omp_get_thread_num();
			dump[i] = graph.dijkstra_independent_to_all_threaded(v[i].id-1,destinations, thread);//return traveltime to all unique targets
		}
	}
	end = clock();
	double time = difftime(end, start) / CLOCKS_PER_SEC;
	cout << endl << "writing time independent travel times to output file after: " << time << " seconds" << endl;
	//write to txt
	ofstream output;
	FILE* fp = NULL;
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	string filename = "titt" + to_string(maxvertices) + ".TXT";
	output.open(filepath+filename, ios::out);
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int j = 0; j < maxvertices; ++j)
		{
			output << dump[i][mapper[j]] << ";";//convert j index to mapped target index
		}
		output << "\n";
	}
	output.close();
}//end construct time independent traveltime

void Ins::construct_time_dependent_traveltime(Graph& graph)
{
	vector<vector<vector<double>>> dump;
	cout << "Constructing time-dependent travel time" << endl;
	clock_t start, end;
	start = clock();
	dump.resize(maxvertices);//no push back in parallel zone
	vector<int>mapper(maxvertices, 0);
	vector<int>destinations;
	for (int i = 0; i < maxvertices; ++i)
	{
		if (find(destinations.begin(), destinations.end(), v[i].id - 1) == destinations.end())//unique beindex
		{
			destinations.push_back(v[i].id - 1);
			mapper[i] = int(destinations.size()) - 1;
		}
		else//non-unique beindex
		{
			//map this vertex to the first occurence in the desination vector
			mapper[i] = int(find(destinations.begin(), destinations.end(), v[i].id - 1) - destinations.begin());
		}
	}
	//start parallel session
	#pragma omp parallel num_threads(12)
	{
		#pragma omp for nowait
		for (int i = 0; i < maxvertices; ++i)
		{
			int thread = omp_get_thread_num();
			dump[i].resize(maxtimeslots);
			for (int t = 0; t < maxtimeslots; ++t)
			{
				dump[i][t] = graph.dijkstra_dependent_to_all_threaded(v[i].id-1,destinations,time_periods[t],thread);
			}// for all timeslots
		}//for all vertices
	}//end pragma parallel
	end = clock();
	double time = difftime(end, start) / CLOCKS_PER_SEC;
	cout << endl << "writing time-dependent traveltime to output file after: " << time << " seconds" << endl;
	ofstream output;
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	string filename = "tt" + to_string(maxvertices) + ".TXT";
	output.open(filepath + filename, ios::out);
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int j = 0; j < maxvertices; ++j)
		{
			for (int t = 0; t <maxtimeslots; ++t)
			{
				output << dump[i][t][mapper[j]] << ";";
			}
			output << "\n";
		}
	}
	output.close(); 
}

void Ins::create_neighbourhood(std::string path, std::string name)
{
	const int N = maxvertices - 2;                 // regular vertices
	const int enddepot = maxvertices - 1;

	// Global bounds & shaping
	const int   K_min = 50;
	const int   K_max = 200;                       // a bit higher for large sets
	const double beta = 3.0;
	const double gamma = 2.0;
	const double eps = 1e-9;

	struct Cand { int j; double s; };              // smaller s is better

	auto clampK = [&](int k) { return std::min(K_max, std::max(K_min, k)); };

	auto keep_top = [](std::vector<Cand>& v, int k) 
		{
		if ((int)v.size() > k) {
			std::nth_element(v.begin(), v.begin() + k, v.end(),
				[](const Cand& a, const Cand& b) { return a.s < b.s; });
			v.resize(k);
		}
		std::sort(v.begin(), v.end(),
			[](const Cand& a, const Cand& b) { return a.s < b.s; });
		};

	// Resize depot containers (we’ll fill later)
	v[enddepot].nb.resize(maxtours);
	v[enddepot].nbi.resize(maxtours);

	// ---------- Per-tour loop ----------
	for (int d = 0; d < maxtours; ++d)
	{
		// -------- Phase 1: precompute top-M (B1) successors for each origin u --------
		const int M_MUTUAL = 20;                   // j keeps i if i is in j’s top-M by B1

		// For membership tests we’ll keep sorted vectors (binary_search)
		std::vector<std::vector<int>> topM_out(maxvertices);

		#pragma omp parallel for schedule(guided)
		for (int u = 0; u < maxvertices - 1; ++u)  // skip depot as origin
		{
			std::vector<Cand> b1_all;
			b1_all.reserve(N);

			// Build feasibility pool for u (reachable j with j→v feasible)
			for (int j = 1; j < maxvertices - 1; ++j) {
				if (j == u) continue;

				// earliest arrival at j leaving u at LTW(u)+serv(u)
				double arr_j = v[u].LTW[d] + v[u].serv + v[u].con[j]->determin;
				if (arr_j < v[j].LTW[d]) arr_j = v[j].LTW[d];
				if (arr_j > v[j].UTW[d]) continue;

				// depart j after service and reach depot in time?
				double dep_j = arr_j + v[j].serv;
				if (dep_j + v[j].con[enddepot]->determin > v[enddepot].UTW[d]) continue;

				// Primary metric B1
				const double tt_uj = v[u].con[j]->determin;
				const double tt_jD = v[j].con[enddepot]->determin;
				const double scorej = (double)v[j].score;
				double s1 = (tt_uj + 0.5 * tt_jD + eps) / (scorej + eps);
				b1_all.push_back({ j, s1 });
			}

			keep_top(b1_all, std::min(M_MUTUAL, (int)b1_all.size()));

			// store sorted list of j for membership checks
			auto& out = topM_out[u];
			out.reserve(b1_all.size());
			for (auto& c : b1_all) out.push_back(c.j);
			std::sort(out.begin(), out.end());
		}

		// -------- Phase 2: build neighborhoods with diversified buckets + mutuals --------
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < maxvertices - 1; ++i)  // skip depot as origin
		{
			// Make sure containers exist
			v[i].nb.resize(maxtours);
			v[i].nbi.resize(maxtours);

			// --------- Feasibility pool ---------
			std::vector<int> pool; pool.reserve(N);
			for (int j = 1; j < maxvertices - 1; ++j) 
			{
				if (j == i) continue;

				double arr_j = v[i].LTW[d] + v[i].serv + v[i].con[j]->determin;
				if (arr_j < v[j].LTW[d]) arr_j = v[j].LTW[d];
				if (arr_j > v[j].UTW[d]) continue;

				double dep_j = arr_j + v[j].serv;
				if (dep_j + v[j].con[enddepot]->determin > v[enddepot].UTW[d]) continue;

				pool.push_back(j);
			}

			//start depot has all vertices as neighbours
			if (i == 0)
			{
				v[i].nb[d].clear();
				v[i].nb[d].reserve(pool.size() + 1);
				for (int j : pool) v[i].nb[d].push_back(&v[j]);
				v[i].nb[d].push_back(&v[enddepot]); // always include end depot

				v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
				v[i].nbi[d].reset();
				for (auto* pj : v[i].nb[d]) v[i].nbi[d][pj->index] = true;
				v[i].nbi[d][i] = true; // self for swap
				continue; // skip trimming/bucketing/mutuals for depot
			}
			
			const int Nreach = (int)pool.size();
			int K_target = clampK((int)std::ceil(beta * std::sqrt(std::max(1, Nreach)) + gamma * std::log1p(Nreach)));
			if (Nreach == 0) 
			{
				v[i].nb[d].clear();
				v[i].nb[d].push_back(&v[enddepot]);
				v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
				v[i].nbi[d].reset();
				v[i].nbi[d][enddepot] = true;
				v[i].nbi[d][i] = true;
				continue;
			}

			// --------- Buckets (B1..B5) as before ---------
			std::vector<Cand> b1, b2, b3, b4, b5;
			b1.reserve(std::min(Nreach, K_target * 3));
			b2.reserve(std::min(Nreach, K_target * 2));
			b3.reserve(std::min(Nreach, K_target * 2));
			b4.reserve(std::min(Nreach, K_target * 2));
			b5.reserve(std::min(Nreach, K_target));

			const double LTWi = v[i].LTW[d];
			const double servi = v[i].serv;

			for (int j : pool) {
				const double tt_ij = v[i].con[j]->determin;
				const double tt_jD = v[j].con[enddepot]->determin;
				const double scorej = (double)v[j].score;

				// B1
				double s1 = (tt_ij + 0.5 * tt_jD + eps) / (scorej + eps);
				b1.push_back({ j, s1 });

				// B2
				double s2 = (tt_ij + eps) / (scorej + eps);
				b2.push_back({ j, s2 });

				// B3 (time-window slack): store negative so larger slack first
				double arr_j = LTWi + servi + tt_ij;
				if (arr_j < v[j].LTW[d]) arr_j = v[j].LTW[d];
				double slack = std::max(0.0, v[j].UTW[d] - arr_j);
				b3.push_back({ j, -slack });
			}

			auto cap = [&](int portion) { return std::min(portion, (int)pool.size()); };
			int k1 = cap((int)std::round(0.90 * K_target));
			int k2 = cap((int)std::round(0.08 * K_target));
			int k3 = cap((int)std::round(0.02 * K_target));

			keep_top(b1, k1); keep_top(b2, k2); keep_top(b3, k3);

			// --------- Union & dedupe of buckets ---------
			std::vector<Cand> uni;
			uni.reserve(b1.size() + b2.size() + b3.size() + b4.size() + b5.size());
			auto add_all = [&](const std::vector<Cand>& src) { uni.insert(uni.end(), src.begin(), src.end()); };
			add_all(b1); add_all(b2); add_all(b3); add_all(b4); add_all(b5);

			// --------- Mutual neighbor additions ---------
			// Add any j where i ∈ topM_out[j]
			for (int j : pool) {
				const auto& list = topM_out[j];
				if (!list.empty() && std::binary_search(list.begin(), list.end(), i)) {
					// Push with B1 score so later trimming uses the same metric
					const double tt_ij = v[i].con[j]->determin;
					const double tt_jD = v[j].con[enddepot]->determin;
					const double scorej = (double)v[j].score;
					double s1 = (tt_ij + 0.5 * tt_jD + eps) / (scorej + eps);
					uni.push_back({ j, s1 });
				}
			}

			// Consolidate: best score per j
			std::sort(uni.begin(), uni.end(), [](const Cand& a, const Cand& b) {
				if (a.j != b.j) return a.j < b.j;
				return a.s < b.s;
				});
			int w = 0;
			for (int r = 0; r < (int)uni.size(); ) {
				int jj = uni[r].j;
				double best = uni[r].s;
				int r2 = r + 1;
				while (r2 < (int)uni.size() && uni[r2].j == jj) {
					if (uni[r2].s < best) best = uni[r2].s;
					++r2;
				}
				uni[w++] = { jj, best };
				r = r2;
			}
			uni.resize(w);

			// Final trimming by B1 to K_target
			keep_top(uni, std::min(K_target, (int)uni.size()));

			// --------- Materialize nb + nbi ---------
			v[i].nb[d].clear();
			v[i].nb[d].reserve(uni.size() + 1);
			for (const auto& c : uni) v[i].nb[d].push_back(&v[c.j]);
			v[i].nb[d].push_back(&v[enddepot]); // always include end depot

			v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
			v[i].nbi[d].reset();
			for (auto* pj : v[i].nb[d]) v[i].nbi[d][pj->index] = true;
			v[i].nbi[d][i] = true; // self for swap
		}

		//end depot neighborhood for this tour
		v[enddepot].nb[d].clear();
		v[enddepot].nb[d].push_back(&v[enddepot]);
		v[enddepot].nbi[d].resize(maxvertices);  
		v[enddepot].nbi[d].reset();
		v[enddepot].nbi[d].set(enddepot);           
	}

	// Write to file
	std::ofstream file(path + "nb" + name);
	for (int i = 0; i < maxvertices; ++i) {
		for (int d = 0; d < maxtours; ++d) {
			file << v[i].nb[d].size() << "\n";
			for (int j = 0; j < (int)v[i].nb[d].size(); ++j)
				file << v[i].nb[d][j]->index << ";";
			file << "\n";
		}
	}
	file.close();
	
	/*
	//full neighborhood (for testing)
	for (int i = 0; i < maxvertices - 1; ++i)  // skip depot as origin
	{
		// Make sure containers exist
		v[i].nb.resize(maxtours);
		v[i].nbi.resize(maxtours);
		for (int d = 0; d < maxtours; ++d)
		{
			v[i].nb[d].clear();
			v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
			v[i].nbi[d].reset();
			for (int j = 0; j < maxvertices; ++j)
			{
				if (i != j)
				{
					v[i].nb[d].push_back(&v[j]);
				}
				v[i].nbi[d][j] = true;
			}
		}
	}
	*/
}

void Ins::create_neighbourhood_simple(std::string path, std::string name)
{
	const int enddepot = maxvertices - 1;
	const int K = 50;              // fixed neighborhood size
	const double eps = 1e-9;

	struct Cand
	{
		int j;
		double s; // smaller is better
	};

	// make sure depot containers exist
	v[enddepot].nb.resize(maxtours);
	v[enddepot].nbi.resize(maxtours);

	for (int d = 0; d < maxtours; ++d)
	{
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < maxvertices - 1; ++i) // skip end depot as origin
		{
			v[i].nb.resize(maxtours);
			v[i].nbi.resize(maxtours);

			std::vector<Cand> cand;
			cand.reserve(maxvertices);

			// ---------- Build feasible pool ----------
			for (int j = 1; j < maxvertices - 1; ++j)
			{
				if (j == i) continue;

				double arr_j = v[i].LTW[d] + v[i].serv + v[i].con[j]->determin;
				if (arr_j < v[j].LTW[d]) arr_j = v[j].LTW[d];
				if (arr_j > v[j].UTW[d]) continue;

				double dep_j = arr_j + v[j].serv;
				if (dep_j + v[j].con[enddepot]->determin > v[enddepot].UTW[d]) continue;

				// simple score: prize per travel with depot look-ahead
				const double tt_ij = v[i].con[j]->determin;
				const double tt_jD = v[j].con[enddepot]->determin;
				const double scorej = (double)v[j].score;

				double s = (tt_ij + 0.5 * tt_jD + eps) / (scorej + eps);
				cand.push_back({ j, s });
			}

			// ---------- Start depot: keep all feasible ----------
			if (i == 0)
			{
				std::sort(cand.begin(), cand.end(),
					[](const Cand& a, const Cand& b) { return a.s < b.s; });

				v[i].nb[d].clear();
				v[i].nb[d].reserve(cand.size() + 1);

				for (const auto& c : cand)
					v[i].nb[d].push_back(&v[c.j]);

				v[i].nb[d].push_back(&v[enddepot]);

				v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
				v[i].nbi[d].reset();
				for (auto* pj : v[i].nb[d]) v[i].nbi[d][pj->index] = true;
				v[i].nbi[d][i] = true;

				continue;
			}

			// ---------- Other vertices: keep top K ----------
			if ((int)cand.size() > K)
			{
				std::nth_element(cand.begin(), cand.begin() + K, cand.end(),
					[](const Cand& a, const Cand& b) { return a.s < b.s; });
				cand.resize(K);
			}

			std::sort(cand.begin(), cand.end(),
				[](const Cand& a, const Cand& b) { return a.s < b.s; });

			v[i].nb[d].clear();
			v[i].nb[d].reserve(cand.size() + 1);

			for (const auto& c : cand)
				v[i].nb[d].push_back(&v[c.j]);

			v[i].nb[d].push_back(&v[enddepot]);

			v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
			v[i].nbi[d].reset();
			for (auto* pj : v[i].nb[d]) v[i].nbi[d][pj->index] = true;
			v[i].nbi[d][i] = true; // self for swap
		}

		// ---------- End depot ----------
		v[enddepot].nb[d].clear();
		v[enddepot].nb[d].push_back(&v[enddepot]);

		v[enddepot].nbi[d].resize(maxvertices);
		v[enddepot].nbi[d].reset();
		v[enddepot].nbi[d].set(enddepot);
	}

	// ---------- Write to file ----------
	std::ofstream file(path + "nb" + name);
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int d = 0; d < maxtours; ++d)
		{
			file << v[i].nb[d].size() << "\n";
			for (int j = 0; j < (int)v[i].nb[d].size(); ++j)
				file << v[i].nb[d][j]->index << ";";
			file << "\n";
		}
	}
	file.close();
}


void Ins::read_neighbourhood(string path,string name)
{
	//read nearest bemobile node
	ifstream file;
	file.open(path+"nb"+name);
	if (file.is_open())
	{
		string value;
		bool stop = false;
		for (int i = 0; i < maxvertices; ++i)//for all regular vertices
		{
			v[i].nbi.resize(maxtours);
			v[i].nb.resize(maxtours);
			for (int d = 0; d < maxtours; ++d)
			{
				getline(file, value, '\n');
				int size = stoi(value);//read amount of neighbours
				v[i].nb[d].resize(size);
				v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
				v[i].nbi[d].set(0);//set bitset to zero for all vertices
				for (int j = 0; j < size; ++j)
				{
					getline(file, value, ';');
					int index = stoi(value);
					v[i].nb[d][j] = &v[index];
					v[i].nbi[d][index] = true;
				}
				getline(file, value, '\n');
				v[i].nbi[d][i] = true;
			}
		}
		file.close();
	}
	else
	{
		cout << "can not find neighborhood file" << endl;
		create_neighbourhood(path,name);
	}
}

void Ins::alter_instance()
{
	breakdur = 0.0;
	
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int j = 0; j < maxvertices; ++j)
		{
			for (int t = 0; t < maxtimeslots; ++t)
			{
				if (t == maxtimeslots - 1)
				{
					v[i].con[j]->mu[t]=0.0;
					v[i].con[j]->oneplusmu[t] = 1.0;
					v[i].con[j]->nu[t]=v[i].con[j]->determin;
				}
				else
				{
					v[i].con[j]->mu[t]=0.0;
					v[i].con[j]->oneplusmu[t] = 1.0;
					v[i].con[j]->nu[t]=v[i].con[j]->determin;
				}
			}
		}
	}
	
}

void Ins::unalter_instance()
{
	breakdur = 0.75;
	
	read_time_dependent_traveltime();
	
}
