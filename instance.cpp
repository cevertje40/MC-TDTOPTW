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
		cout << "reading time-independent travel time" << endl;
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
	cout << "reading time-dependent travel time" << endl;
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
			auto& nu = v[i].con[j]->nu;
			mu.resize(maxtimeslots);
			nu.resize(maxtimeslots);
			//double mintraveltime = DBL_MAX;
			for (int t = 0; t < maxtimeslots; ++t) 
			{
				//mintraveltime = min(mintraveltime, dump[t]);
				if (t == maxtimeslots - 1) {
					mu[t] = 0.0;
					nu[t] = dump[t];
				}
				else {
					double delta_time = time_periods[t + 1] - time_periods[t];
					mu[t] = (dump[t + 1] - dump[t]) / delta_time;
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
				//calculate euclidean distance
				c[counter].determin = sqrt(pow(x_coordinates[i] - x_coordinates[j], 2) + pow(y_coordinates[i] - y_coordinates[j], 2));
				c[counter].mu.resize(maxtimeslots);
				c[counter].nu.resize(maxtimeslots);
				for (int t = 0; t < maxtimeslots; ++t)
				{
					c[counter].mu[t] = 0.0;
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


enum class NbScoreKey {
	TimePerReward,              // A
	DepotTimePerReward,         // B
	SlackAware,                 // C
	DeadlineBias,               // D
	Composite                   // E
};

struct NbScoreParams {
	double alpha = 0.6; // composite blend
	double phi = 0.25; // slack penalty weight
	double psi = 0.2;  // deadline penalty weight
	double SLmin = 0.0001; // slack scale (hours)
};

inline double nb_score(
	const Ins::Vertex& vi, const Ins::Vertex& vj,
	int d, int endIdx, double LATd,
	NbScoreKey key, const NbScoreParams& P)
{
	const double det_ij = vi.con[vj.index]->determin;
	const double det_jv = vj.con[endIdx]->determin;
	const double serv_j = vj.serv;
	const double score_j = max(1.0, (double)vj.score);

	const double T_ij = det_ij + serv_j;
	const double T_ijv = det_ij + serv_j + det_jv;

	const double slack = max(0.0, vj.UTW[d] - vj.LTW[d]);
	const double slack_factor = 1.0 + P.phi * (1.0 - min(1.0, slack / P.SLmin));
	const double urgency = (LATd > 0.0) ? max(0.0, (LATd - vj.UTW[d])) / LATd : 0.0;

	switch (key) {
	case NbScoreKey::TimePerReward:        // A
		return T_ij / score_j;

	case NbScoreKey::DepotTimePerReward:   // B
		return T_ijv / score_j;

	case NbScoreKey::SlackAware:           // C
		return (T_ij / score_j) * slack_factor;

	case NbScoreKey::DeadlineBias:         // D
		return (T_ij / score_j) * (1.0 + P.psi * urgency);

	case NbScoreKey::Composite:            // E
	default: {
		const double blend = P.alpha * T_ij + (1.0 - P.alpha) * T_ijv;
		return (blend / score_j) * slack_factor;
	}
	}
}



void Ins::create_neighbourhood(string path, string name, int amnt_nb)
{
	#pragma omp parallel
	{//start parallel session
		#pragma omp for nowait
		for (int i = 0; i < maxvertices - 1; ++i)//for all regular vertices
		{
			v[i].nb.resize(maxtours);
			v[i].nbi.resize(maxtours);
			for (int d = 0; d < maxtours; ++d)
			{
				//calculate neighbourhood potential
				vector<double> score;
				for (int j = 1; j < maxvertices - 1; ++j)
				{
					if (i != j)
					{
						if (v[i].LTW[d] + v[i].serv + v[i].con[j]->determin <= v[j].UTW[d])
						{
							double arr_j = v[i].LTW[d] + v[i].serv + v[i].con[j]->determin;
							if (arr_j < v[j].LTW[d]) 
							{
								arr_j = v[j].LTW[d]; // wait for j to open if needed
							}
							// Earliest depart from j after service
							double dep_j = arr_j + v[j].serv;
							double eta_depot = dep_j + v[j].con[maxvertices - 1]->determin;
							//departing at ltw +serv from i you have to be able to reach depot from j before the depots' utw
							if (eta_depot <= v[maxvertices-1].UTW[d])
							{
								score.push_back(max(1.0, v[i].con[j]->determin + v[j].serv) / v[j].score);
								v[i].nb[d].push_back(&v[j]);
								//NbScoreParams P; // defaults OK; tweak if desired
								//const auto key = NbScoreKey::Composite; // choose one
								//const double s = nb_score(v[i], v[j], d, maxvertices - 1, t[d].LAT, key, P);
								//score.push_back(s);
								//v[i].nb[d].push_back(&v[j]);
							}
						}
					}
				}//end for j
				// sort based on score potential
				bool unsorted = true;
				while (unsorted)
				{
					unsorted = false;
					for (int j = 0; j < (int)score.size() - 1; ++j)
					{
						if (score[j] > score[j + 1])
						{
							Vertex* tempnb = v[i].nb[d][j];
							double temps = score[j];
							v[i].nb[d][j] = v[i].nb[d][j + 1];
							score[j] = score[j + 1];
							score[j + 1] = temps;
							v[i].nb[d][j + 1] = tempnb;
							unsorted = true;
						}
					}//end for
				}//end while
				//select top elements
				if ((int)v[i].nb[d].size() > amnt_nb)
				{
					v[i].nb[d].erase(v[i].nb[d].begin() + amnt_nb, v[i].nb[d].end());
				}
				//automatically add the end depot
				v[i].nb[d].push_back(&v[maxvertices - 1]);
				//indexed list aanmaken
				v[i].nbi[d] = boost::dynamic_bitset<>(maxvertices);
				v[i].nbi[d].reset();//sets all bits to false
				for (int x = 0; x < (int)v[i].nb[d].size(); ++x)
				{
					v[i].nbi[d][v[i].nb[d][x]->index] = true;
				}
				//for a swap operation each vertex is neighbour of itself
				v[i].nbi[d][i] = true;
			}//end for d
		}//end for i
	}//end parallel
	//and enddepot to enddepot
	v[maxvertices - 1].nb.resize(maxtours);
	v[maxvertices - 1].nbi.resize(maxtours);
	for (int d = 0; d < maxtours; ++d)
	{
		v[maxvertices - 1].nbi[d] = boost::dynamic_bitset<>(maxvertices);
		//v[maxvertices - 1].nbi[d].set(0);
		v[maxvertices - 1].nb[d].push_back(&v[maxvertices - 1]);
		v[maxvertices - 1].nbi[d][maxvertices - 1] = true;
	}
	ofstream file;
	file.open(path+"nb" + name);
	for (int i = 0; i < maxvertices; ++i)//for all regular vertices
	{
		for (int d = 0; d < maxtours; ++d)
		{
			file << v[i].nb[d].size() << "\n";
			for (int j = 0; j < v[i].nb[d].size(); ++j)
			{
				file << v[i].nb[d][j]->index << ";";
			}
			file << "\n";
		}
	}
	file.close();
}//end neighbourhood

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
		create_neighbourhood(path,name,50);
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


inline int Ins::find_t(double time)
{
	int t = (int)floor((time - time_periods[0]) / 0.25);//when you change the time unit this has to change too
	return min(55,t);
}

double Ins::travel_time(Connec* c, double start)
{
	int t = find_t(start);
	double traveltime = c->nu[t] + start * c->mu[t];
	return traveltime;
}

double Ins::arrival_time(Connec* c, double start)
{
	int t = find_t(start);
	double arrivaltime = c->nu[t] + (start)*c->mu[t] + start;
	return arrivaltime;
}

double Ins::departure_time(Connec* c, double arrivaltime)
{
	int t = find_t(arrivaltime);
	double departuretime = (arrivaltime - c->nu[t]) / (1 + c->mu[t]);
	while ((time_periods[t] > departuretime) || (departuretime > time_periods[t + 1]))
	{
		--t;
		departuretime = (arrivaltime - c->nu[t]) / (1 + c->mu[t]);
	}
	return departuretime;
}
