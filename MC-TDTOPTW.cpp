// MC-TDTOPTW.cpp

#include "solver.h"

using namespace std;

void create_dataset()
{
	//user input
	vector<vector<int>>vertexid(3, vector<int>());
	vertexid[0] = {42,18468,26501,15725,29359,24465,28146,16828,492,11943,5437,14605,154,12383,18717,19896,21727,11539,19913,26300};
	vertexid[1] = {11174,10467,21660,26440,20025,29511,20650,8314,28023,14019,9906,7392,3626,4415,25825,25875,20160,28071,28298,8178,32271,2669,13986,8481,7628,4100,2626,1925,29973,14182,27433,27594,13032,143,31287,7901,8361,30975,29171,30834,25761,4668,12551,13695,21625,2126,21695,26303,22467,22594 };
	vertexid[2] = {27466,20268,19794,25473,22831,28443,13878,703,1382,8824,8024,16596,2328,31311,11059,9488,32529,2259,9861,21287,8611,7129,5842,3504,24866,1882,22751,18599,2662,32757,20279,19436,32076,1387,8361,26049,29493,23841,1736,11600,21893,7329,11370,21795,9253,17433,7209,3498,27650,26842,16101,30649,19852,28634,27201,9991,4920,22579,32545,13488,22526,5539,6194,25012,15835,31498,18530,18806,13393,13550,26980,9278,20194,21498,31277,6583,11160,26490,3450,9073,27009,10209,18504,32608,12075,12612,28762,12891,16684,19933,2742,6814,10397,20616,2600,4681,27033,32585,3518,8671};
	vector<double> tmaxarray{7.0,8.0,9.0};
	vector<double> twseverity{0.8,0.6,0.4};
	vector<int> tours{2,3,4};
	double t_zero = 6;
	double servmean = 0.33;
	double servsd = 0.1;
	double servkappa=pow(servmean, 2) / pow(servsd, 2);
	double servtheta= pow(servsd, 2) / servmean;
	double weightmean = 3.0;
	double weightsd = 2.0;
	double weightkappa = pow(weightmean, 2) / pow(weightsd, 2);
	double weighttheta = pow(weightsd, 2) / weightmean;
	double volmean = 5.0;
	double volsd = 2.0;
	double volkappa = pow(volmean, 2) / pow(volsd, 2);
	double voltheta = pow(volsd, 2) / volmean;
	int max_score = 40;
	double breakstart = t_zero + 3;
	double breakend = t_zero + 7;
	double breakdur = 0.75;
	//automatic code
	default_random_engine generator;
	generator.seed((unsigned)time(0));
	gamma_distribution<double> servicedist(servkappa,servtheta);
	gamma_distribution<double> voldist(volkappa, voltheta);
	gamma_distribution<double> weightdist(weightkappa,weighttheta);
	cout << servicedist(generator) << endl;
	for (int v = 0; v < vertexid.size(); ++v)//for all 3 maxvertices
	{
		for (int t = 0; t < tours.size(); ++t)
		{
			for (int tm = 0; tm < tmaxarray.size(); ++tm)
			{
				for (int tw = 0; tw < twseverity.size(); ++tw)
				{
					int maxvertices = int(vertexid[v].size());
					vector<int> ids(maxvertices);
					vector<int> scores(maxvertices);
					vector<double>services(maxvertices);
					vector<double>volumes(maxvertices);
					vector<double>weights(maxvertices);
					int maxtours = tours[t];
					vector<vector<double>>ltws(maxvertices,vector<double>(maxtours,0.0));
					vector<vector<double>>utws(maxvertices,vector<double>(maxtours,0.0));
					double T_max = tmaxarray[tm];
					double W_max = 22;//tons standard truck
					double V_max = 47;//m3 standard truck
					for (int i = 0; i < maxvertices; ++i)
					{
						ids[i] = vertexid[v][i];
						if ((i == 0) || (i == maxvertices - 1))//start & end vertex
						{
							scores[i] = 0;
							services[i] =0;
							volumes[i] = 0;
							weights[i] = 0;
							for (int b = 0; b < maxtours; ++b)
							{
								ltws[i][b] = t_zero;
								utws[i][b] = t_zero + T_max;
							}
						}
						else
						{//regular vertex
							scores[i] = 1 + rand() % (max_score - 1);
							services[i] = servicedist(generator);
							weights[i] = 0.5 + weightdist(generator);
							volumes[i] = 2+voldist(generator);
							int severity = int(twseverity[tw] * T_max);
							for (int b = 0; b < maxtours; ++b)
							{
								int dividora = int(T_max - severity);
								ltws[i][b] = t_zero + rand() % dividora;
								int dividorb = int(t_zero + T_max - (severity + ltws[i][b]) + 1);
								utws[i][b] = ltws[i][b] + severity + rand() % dividorb;
							}
						}
					}
					//write to file
					ofstream output;
					cout.precision(3);
					string name= to_string(maxvertices)+ "." + to_string(t + 1) + "." + to_string(tm + 1) + "." + to_string(tw + 1) +".txt";
					output.open(name, ios::out);
					output << maxvertices << '\n';
					output << maxtours << '\n';
					output << breakdur << '\t' << breakstart << '\t'<< breakend << '\n';
					for (int a = 0; a < vertexid[v].size(); ++a)
					{
						output << ids[a] << '\t' << scores[a] << setprecision(2) << '\t' << services[a] << '\t' << weights[a] << '\t' << volumes[a] << '\t';
						for (int b = 0; b < maxtours; ++b)
						{
							output << ltws[a][b] << '\t' << utws[a][b] << '\t';
						}
						output << '\n';
					}
					for (int b = 0; b < maxtours; ++b)
					{
						output << T_max << '\t' << W_max << '\t' << V_max << '\n';
					}
					output.close();
					
				}//for all tw values
			}//for all tmax values
		}//for all tour values
	}//for all maxvertex size values
}

list <string> read_dataset(string filename)
{//reads in all the dataset names
	list <string> files; //filenames list with iterator
	ifstream ifs;
	ifs.open(filename, ifstream::in);
	if (ifs.is_open())
	{
		string line;
		while (getline(ifs, line))//read 1 full line
		{
			stringstream str(line);//store line as stringstream
			files.push_back(line);
		}
		ifs.close();
	}
	else
	{
		printf("\ninput error in filenames file");
	}
	return files;
}

void solve_dataset(int testruns)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(std::cin, filename);
	if (filename.size() == 0)
	{
		filename = "all.txt";
	}
	list<string> files = read_dataset(filename);
	vector<Res> resdataset;
	list <string>::iterator it;
	for (it = files.begin(); it != files.end(); ++it)
	{
		cout << *it << endl;
		Ins instance(*it);
		instance.read_neighbourhood();
		Aco acs(instance, 1, 3, 0.01, 20, 10000, 0.25, 0.05);
		resdataset.push_back(acs.solve());
	}
}

void debug_instance()
{
	vector<Res> resdataset;
	Ins instance("20.1.3.2.txt");
	instance.read_neighbourhood();
	Aco acs(instance, 1, 3, 0.01, 20, 10000, 0.25, 0.05);
	resdataset.push_back(acs.solve());
}

int main()
{
	solve_dataset(1);
	//debug_instance();
}
