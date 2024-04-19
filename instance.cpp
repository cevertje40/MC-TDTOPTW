#include "instance.h"

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
		cout << red << "could not open time-independent travel time file" << endl;
	}

}

void Ins::read_time_dependent_traveltime()
{
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	ifstream tt;
	tt.open(filepath + "tt" + to_string(maxvertices) + ".TXT", ifstream::in);
	if (tt.is_open())
	{
		cout << "reading time-dependent travel time" << endl;
		for (int i = 0; i < maxvertices; ++i)
		{
			for (int j = 0; j < maxvertices; ++j)
			{
				vector<double> dump(maxtimeslots);
				string line;
				getline(tt, line);
				stringstream str(line);
				for (int t = 0; t < maxtimeslots; ++t)
				{
					str >> dump[t];
					str.ignore();
				}
				for (int t = 0; t < maxtimeslots; ++t)
				{
					if (t == maxtimeslots - 1)
					{
						v[i].con[j]->mu.push_back(double(dump[t] - dump[t]) / (time_periods[t + 1] - time_periods[t]));
						v[i].con[j]->nu.push_back(double((dump[t]) - double(v[i].con[j]->mu[t] * time_periods[t])));
					}
					else
					{
						v[i].con[j]->mu.push_back(double(dump[t + 1] - dump[t]) / (time_periods[t + 1] - time_periods[t]));
						v[i].con[j]->nu.push_back(double((dump[t]) - double(v[i].con[j]->mu[t] * time_periods[t])));
					}
				}
			}
		}
		tt.close();
	}
	else
	{
		cout << red << "could not open time-dependent travel time file" << endl;
	}
}

Ins::Ins(string filename)
{
	name = filename;
	//read in vertex and tour information from txt file and populate v and t objects
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	FILE* file = NULL;
	ifstream ifs;
	ifs.open(filepath + filename, ifstream::in);
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
	//read in travel time information
	c.clear();
	c.resize(maxvertices * maxvertices);
	read_time_independent_traveltime();
	read_time_dependent_traveltime();

}

void Ins::construct_time_independent_traveltime(Graph& graph)
{
	clock_t start, end;
	start = clock();
	vector<int> targets;
	vector<vector<double>> dump;
	for (int i = 0; i < maxvertices; ++i)
	{
		targets.push_back(v[i].id-1);//todo only push back unique beindex
		dump.push_back(vector<double>());
	}
	cout << "construct time independent travel time" << endl;
	#pragma omp parallel num_threads(12)
	{//start parallel session
		#pragma omp for nowait
		for (int i = 0; i < maxvertices; ++i)
		{
			int thread = omp_get_thread_num();
			dump[i] = graph.dijkstra_independent_to_all_threaded(v[i].id-1, targets, thread);
		}
	}
	end = clock();
	double time = difftime(end, start) / CLOCKS_PER_SEC;

	//write to txt
	FILE* fp = NULL;
	char filepath[125] = "..\\..\\datasets\\MCTDTOPTW\\";
	char filename[15];
	sprintf_s(filename, sizeof(filename), "titt%d.TXT", maxvertices);
	strcat_s(filepath, filename);
	fopen_s(&fp, filepath, "w");   // open for writing 
	if (fp != NULL)
	{
		cout << endl << "writing time independent travel times to output file after: " << time << " seconds" << endl;
		for (int i = 0; i < maxvertices; ++i)
		{
			for (int j = 0; j < maxvertices; ++j)
			{
				fprintf(fp, "%lf;", dump[i][j]);//output in miliseconds
			}
			fprintf(fp, "\n");
		}
		fclose(fp);    // close the file before ending program 
	}
	else
	{
		cout << red << "error writing time-independent travel times to output file" << endl;
	}
}//end construct time independent traveltime

void Ins::construct_time_dependent_traveltime(Graph& graph)
{
	vector<vector<vector<double>>> dump;
	cout << "Constructing time-dependent travel time" << endl;
	clock_t start, end;
	start = clock();
	//niet in parallele zone push backen
	dump.resize(maxvertices);
	//start parallel session
	#pragma omp parallel num_threads(12)
	{
	#pragma omp for nowait
		for (int i = 0; i < maxvertices; ++i)
		{
			int thread = omp_get_thread_num();
			//obtain neighbours
			vector<int> destinations;
			for (int j = 0; j < maxvertices; ++j)
			{
				destinations.push_back(v[j].id-1);//non neighbours will get infinity
			}
			//define feasible departure time zone
			dump[i].resize(maxtimeslots);
			for (int t = 0; t < maxtimeslots; ++t)
			{
				dump[i][t] = graph.dijkstra_dependent_to_all_threaded(v[i].id-1, destinations,time_periods[t],thread);//j
			}// for all timeslots
		}//for all vertices
	}//end pragma parallel
	end = clock();
	double time = difftime(end, start) / CLOCKS_PER_SEC;

	FILE* file = NULL;
	char filepath[125] = "..\\..\\datasets\\MCTDTOPTW\\";
	char storagename[50];
	sprintf_s(storagename, sizeof(storagename), "tt%d.TXT", maxvertices);
	strcat_s(filepath, storagename);
	fopen_s(&file, filepath, "w");   // open for writing
	if (file != NULL)
	{
		cout << endl << "writing time-dependent traveltime to output file after: " << time << " seconds" << endl;
		for (int i = 0; i < maxvertices; ++i)
		{
			for (int j = 0; j < maxvertices; ++j)
			{
				for (int t = 0; t <maxtimeslots; ++t)
				{
					fprintf(file, "%lf;", dump[i][t][j]);
				}
				fprintf(file, "\n");
			}
		}
		fclose(file);    // close the file before ending program 
	}
	else
	{
		cout << red << "error writing time-dependent traveltime to output file" << endl;
	}
}

void Ins::create_neighbourhood(int amnt_nb)
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
							score.push_back((v[i].con[j]->determin + v[j].serv + v[j].weight + v[j].volume) / v[j].score);
							v[i].nb[d].push_back(&v[j]);
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
				v[i].nbi[d].set(0);//sets all bits to false
				for (int x = 0; x < (int)v[i].nb[d].size(); ++x)
				{
					v[i].nbi[d][v[i].nb[d][x]->index] = true;
				}
				//for a swap operation each vertex is neighbour of itself
				v[i].nbi[d][i] = true;
			}//end for d
		}//end for i
		//and enddepot to enddepot
		v[maxvertices - 1].nb.resize(maxtours);
		v[maxvertices - 1].nbi.resize(maxtours);
		for (int d = 0; d < maxtours; ++d)
		{
			v[maxvertices - 1].nbi[d] = boost::dynamic_bitset<>(maxvertices);
			v[maxvertices - 1].nb[d].push_back(&v[maxvertices - 1]);
			v[maxvertices - 1].nbi[d][maxvertices - 1] = true;
		}
	}//end parallel
	ofstream file;
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	file.open(filepath+"nb" + name);
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

void Ins::read_neighbourhood()
{
	//read nearest bemobile node
	ifstream file;
	string filepath = "..\\..\\datasets\\MCTDTOPTW\\";
	file.open(filepath+"nb"+name);
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
		cout << "error reading neighborhood" << endl;
		create_neighbourhood(45);
	}
}

inline int Ins::find_t(double time)
{
	int t = (int)floor((time - time_periods[0]) / 0.25);//when you change the time unit this has to change too
	return t;
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
