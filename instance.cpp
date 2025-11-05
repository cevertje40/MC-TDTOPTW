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


/*
void Ins::read_time_dependent_traveltime()//with fifo checker
{
	const std::string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	std::ifstream tt(filepath + "tt" + std::to_string(maxvertices) + ".TXT");
	if (!tt.is_open()) {
		std::cout << term::fg(term::Color::red)
			<< "could not open time-dependent travel time file\n";
		return;
	}

	// Tolerances for floating-point checks (hours)
	const double TOL_BOUNDARY = 1e-9;   // boundary FIFO (arr[k+1] + tol >= arr[k])
	const double TOL_SLOPE = 1e-12;  // slope checks (1+mu >= -tol), |mu| <= 1 + tol

	// Stats
	std::uint64_t parse_warn = 0, neg_warn = 0, naninf_warn = 0;
	std::uint64_t fifo_boundary_viol = 0, fifo_withinslot_viol = 0, abs_mu_warn = 0;
	std::uint64_t lines_read = 0;

	for (int i = 0; i < maxvertices; ++i) {
		for (int j = 0; j < maxvertices; ++j) {

			std::string line;
			if (!std::getline(tt, line)) {
				std::cout << term::fg(term::Color::red)
					<< "Unexpected EOF at arc (" << i << "," << j << ")\n";
				return;
			}
			++lines_read;

			// Normalize delimiters: turn ';' ',' '\t' into space so stringstream >> works
			for (char& c : line) {
				if (c == ';' || c == ',' || c == '\t') c = ' ';
			}

			std::stringstream ss(line);
			std::vector<double> ttvals;
			ttvals.reserve(maxtimeslots);
			double val;
			while (ss >> val) ttvals.push_back(val);

			if (static_cast<int>(ttvals.size()) != maxtimeslots) {
				++parse_warn;
				std::cout << term::fg(term::Color::yellow)
					<< "Line " << lines_read << " for arc (" << i << "," << j
					<< "): expected " << maxtimeslots << " values, got "
					<< ttvals.size() << "\n";
				// best-effort: pad/truncate to expected length
				ttvals.resize(maxtimeslots, (ttvals.empty() ? 0.0 : ttvals.back()));
			}

			// Basic value checks: finite and >= 0 (times are in hours)
			for (int k = 0; k < maxtimeslots; ++k) {
				if (!std::isfinite(ttvals[k])) {
					++naninf_warn;
					std::cout << term::fg(term::Color::yellow)
						<< "Non-finite value at arc (" << i << "," << j
						<< "), slot " << k << "\n";
					ttvals[k] = 0.0; // sanitize to keep going
				}
				if (ttvals[k] < 0.0) {
					++neg_warn;
					std::cout << term::fg(term::Color::yellow)
						<< "Negative travel time at arc (" << i << "," << j
						<< "), slot " << k << " : " << ttvals[k] << "\n";
					ttvals[k] = 0.0; // clamp
				}
			}

			// Boundary FIFO check: arrival at boundary must be nondecreasing
			// A_k = time_periods[k] + ttvals[k]
			for (int k = 0; k + 1 < maxtimeslots; ++k) {
				const double Ak = time_periods[k] + ttvals[k];
				const double Ak1 = time_periods[k + 1] + ttvals[k + 1];
				if (Ak1 + TOL_BOUNDARY < Ak) {
					++fifo_boundary_viol;
					// Print first few only to avoid flooding
					if (fifo_boundary_viol <= 5) {
						std::cout << term::fg(term::Color::yellow)
							<< "Boundary FIFO violation on arc (" << i << "," << j<< "), slot " << k << " --> " << (k + 1)<< " : Ak=" << Ak << " > Ak+1=" << Ak1 << "\n";
					}
				}
			}

			// Store coefficients
			auto* edge = v[i].con[j];
			if (!edge) {
				std::cout << term::fg(term::Color::red)
					<< "Null edge pointer at (" << i << "," << j << ")\n";
				return;
			}
			auto& mu = edge->mu;
			auto& oneplusmu = edge->oneplusmu;
			auto& nu = edge->nu;
			mu.resize(maxtimeslots);
			oneplusmu.resize(maxtimeslots);
			nu.resize(maxtimeslots);

			// Compute μ, ν (affine-in-slot model).
			for (int t = 0; t < maxtimeslots; ++t) {
				if (t == maxtimeslots - 1) {
					// last slot: hold travel time constant within the slot
					mu[t] = 0.0;
					oneplusmu[t] = 1.0;
					nu[t] = ttvals[t];
					continue;
				}
				const double dt = time_periods[t + 1] - time_periods[t];
				if (dt <= 0.0) {
					std::cout << term::fg(term::Color::red)
						<< "Non-positive slot width at slot " << t
						<< " (dt=" << dt << ")\n";
					return;
				}

				const double slope = (ttvals[t + 1] - ttvals[t]) / dt;
				mu[t] = slope;
				oneplusmu[t] = 1.0 + slope;
				nu[t] = ttvals[t] - slope * time_periods[t];

				// Within-slot FIFO: 1 + μ >= 0 (allow tiny tolerance)
				if (oneplusmu[t] < -TOL_SLOPE) {
					++fifo_withinslot_viol;
					if (fifo_withinslot_viol <= 5) {
						std::cout << term::fg(term::Color::yellow)
							<< "Within-slot FIFO violated on arc (" << i << "," << j<< "), slot " << t << " : 1+mu=" << oneplusmu[t] << "\n";
					}
				}
				// Empirical check: |μ| <= 1 (your practice)
				if (std::abs(mu[t]) > 1.0 + TOL_SLOPE) {
					++abs_mu_warn;
					if (abs_mu_warn <= 5) {
						std::cout << term::fg(term::Color::yellow)
							<< "Slope magnitude > 1 on arc (" << i << "," << j<< "), slot " << t << " : mu=" << mu[t] << "\n";
					}
				}
			}

			// last slot already filled above (μ=0, ν=ttvals[last])
		}
	}

	tt.close();

	// Summary
	if (parse_warn || naninf_warn || neg_warn || fifo_boundary_viol || fifo_withinslot_viol || abs_mu_warn) {
		std::cout << term::fg(term::Color::yellow)
			<< "[TT checks] lines=" << lines_read
			<< ", parse=" << parse_warn
			<< ", nan/inf=" << naninf_warn
			<< ", negative=" << neg_warn
			<< ", FIFO_boundary=" << fifo_boundary_viol
			<< ", FIFO_within=" << fifo_withinslot_viol
			<< ", |mu|>1=" << abs_mu_warn << "\n";
	}
	else {
		std::cout << term::fg(term::Color::green)
			<< "[TT checks] OK: parsed " << lines_read
			<< " lines; no issues detected.\n";
	}
}
*/

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
	/*
	for (int i = 0; i < maxvertices; ++i)
	{
		for (int j = 0; j < maxvertices; ++j)
		{
			for (int t = 0; t < maxtimeslots; ++t)
			{
				if (t == maxtimeslots - 1)
				{
					v[i].con[j]->mu[t]=0.0;
					v[i].con[j]->nu[t]=v[i].con[j]->determin;
				}
				else
				{
					v[i].con[j]->mu[t]=0.0;
					v[i].con[j]->nu[t]=v[i].con[j]->determin;
				}
			}
		}
	}
	*/
}

void Ins::unalter_instance()
{
	breakdur = 0.75;
	//read_time_dependent_traveltime();
}
