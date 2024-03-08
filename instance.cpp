#include "instance.h"


Ins::Ins(string filename)
{
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
		for (int vertex = 0; vertex < maxvertices; ++vertex)
		{
			getline(ifs, line);//maxvertices lines vertex info
			str = stringstream(line);
			v[vertex].index = vertex;
			str >> v[vertex].id >> v[vertex].score >> v[vertex].serv >> v[vertex].weight >> v[vertex].volume;
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
			t[tour].EDT = 0;
			t[tour].LAT = 0 + t[tour].T_max;
			t[tour].startv = &v[0];
			t[tour].endv = &v[maxvertices - 1];
		}
		ifs.close();
	}
	else
	{
		cout << " can not open MC-TDTOPTW instance file" << endl;
	}
	//read in travel time information from files and store to c objects
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
			dump[i].resize(graph.maxtimeslots);
			for (int t = 0; t < graph.maxtimeslots; ++t)
			{
				dump[i][t] = graph.dijkstra_dependent_to_all_threaded(v[i].id-1, destinations, graph.time_periods[t], thread);//j
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
				for (int t = 0; t < graph.maxtimeslots; ++t)
				{
					fprintf(file, "%lf;", dump[i][t][j]);
				}
			}
		}
		fclose(file);    // close the file before ending program 
	}
	else
	{
		cout << red << "error writing time-dependent traveltime to output file" << endl;
	}
}
