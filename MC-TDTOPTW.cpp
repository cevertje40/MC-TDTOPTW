// MC-TDTOPTW.cpp

#include "solver.h"
using namespace std;

class Instance
{
public:
	string path;
	string filename;
	vector<Res> result;
	int bestknown{0};
	double avgscore{0.0};
	double sdscore{0.0};
	double avggap{0.0};
	double sdgap{0.0};
	double avgtime{0.0 };
	Instance(const string& path, const string& filename, int bestknown) : path(path), filename(filename), bestknown(bestknown){}
	void calculate_statistics()
	{
		avgscore = 0.0;
		sdscore = 0.0;
		avggap = 0.0;
		sdgap = 0.0;
		avgtime = 0.0;
		const int n = (int) result.size();

		if (n == 0) return;//no results yet

		if (n == 1) {// single result: sd = 0
			avgscore = result[0].sol.score;
			avggap = result[0].gap;
			avgtime = result[0].time;
			return;
		}

		for (int i = 0; i < n; ++i)
		{
			avgscore += result[i].sol.score;
			avggap += result[i].gap;
			avgtime += result[i].time;
		}
		avgscore /= n;
		avggap /= n;
		avgtime /= n;
		for (int i = 0; i < n; ++i)
		{
			sdscore += pow((result[i].sol.score - avgscore), 2);
			sdgap += pow((result[i].gap - avggap), 2);
		}
		sdscore = sqrt(sdscore / (n-1));
		sdgap = sqrt(sdgap / (n-1));
	}
};



void create_dataset()
{
	//user input
	vector<vector<int>>vertexid(3, vector<int>());
	vertexid[0] = {42,18468,26501,15725,29359,24465,28146,16828,492,11943,5437,14605,154,12383,18717,19896,21727,11539,19913,26300};
	vertexid[1] = {11174,10467,21660,26440,20025,29511,20650,8314,28023,14019,9906,7392,3626,4415,25825,25875,20160,28071,28298,8178,32271,2669,13986,8481,7628,4100,2626,1925,29973,14182,27433,27594,13032,143,31287,7901,8361,30975,29171,30834,25761,4668,12551,13695,21625,2126,21695,26303,22467,22594 };
	vertexid[2] = {27466,20268,19794,25473,22831,28443,13878,703,1382,8824,8024,16596,2328,31311,11059,9488,32529,2259,9861,21287,8611,7129,5842,3504,24866,1882,22751,18599,2662,32757,20279,19436,32076,1387,8361,26049,29493,23841,1736,11600,21893,7329,11370,21795,9253,17433,7209,3498,27650,26842,16101,30649,19852,28634,27201,9991,4920,22579,32545,13488,22526,5539,6194,25012,15835,31498,18530,18806,13393,13550,26980,9278,20194,21498,31277,6583,11160,26490,3450,9073,27009,10209,18504,32608,12075,12612,28762,12891,16684,19933,2742,6814,10397,20616,2600,4681,27033,32585,3518,8671};
	vector<double> tmaxarray{7.0,8.0,9.0};
	vector<double> twseverity{0.8,0.6,0.4};//lower value means more strict Large, Medium
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
	double W_max = 22;//tons standard truck
	double V_max = 47;//m3 standard truck
	//automatic code
	default_random_engine generator;
	generator.seed((unsigned)time(0));
	gamma_distribution<double> servicedist(servkappa,servtheta);
	gamma_distribution<double> voldist(volkappa, voltheta);
	gamma_distribution<double> weightdist(weightkappa,weighttheta);
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
}//end create dataset

void create_case_dataset()
{
	//user input
	vector<int> vertexid = {41123,313785,13258,385118,3720,364264,157776,166244,387008,386129,386129,160818,327049,402721,44863,270930,7062,19163,372338,373339,18180,244826,22219,56570,362153,5601,212157,16160,371072,326079,371580,416467,42104,174355,175680,765,368836,214376,40173,39510,16424,19052,326308,373972,929,13258,231494,347964,40298,4194,24407,210867,38016,5925,16018,177771,347964,60597,160818,60536,40355,171803,396750,211108,18180,162644,362155,325548,43811,160818,210009,48395,19706,267041,378640,325186,291637,329252,287775,212369,271071,168058,5601,161516,37684,16959,365112,392661,290308,182322,387026,240713,263072,48395,162551,381984,10557,20563,14687,24418,369298,344288,227572,368017,363690,38688,231863,172198,36279,161055,183414,165179,4077,38916,48397,212369,311895,271807,287566,307372,18573,211242,212213,39491,53629,4072,356315,269163,159743,35070,328981,240699,43224,38403,41023,380583,372743,328981,16024,288005,212847,326268,349705,7534,176810,313785,175640,18205,11867,268875,357301,17459,373339,358281,20120,37882,160841,21668,42617,40355,17854,41123};
	vector<double>weightid = {0.0,1.72,0.94,2.72,1.95,0.93,1.4,1.95,0.93,0.31,0.52,1.27,1.44,2.36,0.26,2.25,2.48,0.93,8.85,0.11,1.96,1.61,0.06,1.34,2.87,0.21,1.45,1.56,0.84,1.52,3.84,21.86,1.34,2.97,0.95,1.34,3.24,2.03,0.24,0.77,1.38,0.93,2.8,0.15,2,0.21,3.34,1.85,1.04,2.31,2.49,1,2.62,0.46,6.3,1.44,1.85,0.21,0.11,0.93,0.16,2.43,2.86,1.66,0.94,3.05,7.69,1.34,2.9,0.06,6,0.24,1.64,1.24,0.77,20.95,1.64,3.24,4.5,1.06,4.03,2.48,0.52,3.4,0.61,0.11,2.11,3.34,3.36,2.49,0.69,0.93,5.69,3.26,0.31,0.91,2.63,3.91,0.72,1.02,2.37,35.96,0.56,3.44,0.46,14.03,1.56,0.42,1.2,0.33,28.03,0.93,3.41,2.35,0.93,1.34,0.61,0.43,5.22,0.93,2.34,1.28,2,7.53,0.76,0.65,1.34,2.52,0.93,19.66,0.21,0.93,1.34,0.51,1.34,2.46,1.29,3.49,4.42,1.54,18.23,2.07,8.06,0.79,4.74,0.21,0.8,2.15,4.29,0.62,0.92,6.6,0.51,2.64,3.29,3.11,0.93,0.78,0.93,3.81,0.3,0.0};
	int max_score = 40;
	vector<int> scoreid(vertexid.size(),0);
	int sumscore = 0;
	for (int i = 1; i < vertexid.size()-1; ++i)
	{
		scoreid[i]= 1 + rand() % (max_score - 1);
		sumscore += scoreid[i];
	}
	vector<double> tmaxarray{6.0,8.0};
	vector<double> twseverity{0.8,0.6};//lower value means more strict (Large, Medium)
	vector<int> tours{8,10};
	double t_zero = 6;
	double servmean = 0.2;//deterministic in the case
	double volmean = 1.0;
	double breakstart = t_zero + 3.5;
	double breakend = t_zero + 6;
	double breakdur = 0.75;
	double W_max = 60;//tons standard truck
	double V_max = 20;//m3 standard truck
	//automatic code
	default_random_engine generator;
	generator.seed((unsigned)time(0));
		for (int t = 0; t < tours.size(); ++t)
		{
			for (int tm = 0; tm < tmaxarray.size(); ++tm)
			{
				for (int tw = 0; tw < twseverity.size(); ++tw)
				{
					int maxvertices = int(vertexid.size());
					vector<int> ids(maxvertices);
					vector<int> scores(maxvertices);
					vector<double>services(maxvertices);
					vector<double>volumes(maxvertices);
					vector<double>weights(maxvertices);
					int maxtours = tours[t];
					vector<vector<double>>ltws(maxvertices, vector<double>(maxtours, 0.0));
					vector<vector<double>>utws(maxvertices, vector<double>(maxtours, 0.0));
					double T_max = tmaxarray[tm];
					for (int i = 0; i < maxvertices; ++i)
					{
						ids[i] = vertexid[i];
						if ((i == 0) || (i == maxvertices - 1))//start & end vertex
						{
							scores[i] = 0;
							services[i] = 0;
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
							scores[i] = scoreid[i];
							services[i] = servmean;
							weights[i] = weightid[i];
							volumes[i] = volmean;
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
					string name = to_string(maxvertices) + "." + to_string(t + 1) + "." + to_string(tm + 1) + "." + to_string(tw + 1) + ".txt";
					output.open(name, ios::out);
					output << maxvertices << '\n';
					output << maxtours << '\n';
					output << breakdur << '\t' << breakstart << '\t' << breakend << '\n';
					for (int a = 0; a < vertexid.size(); ++a)
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
}//end create case

vector<Instance> read_dataset(string filename)
{//reads in all the dataset names
	vector<Instance> dataset;
	string path;
	ifstream ifs;
	ifs.open(filename, ifstream::in);
	if (ifs.is_open())
	{
		string line;
		getline(ifs, line);
		stringstream str(line);
		path = line;
		while (getline(ifs, line))//read 1 full line
		{
			str=stringstream(line);//store line as stringstream
			string name;
			int bestknown;
			str>> name;
			str>> bestknown;
			dataset.push_back(Instance(path, name, bestknown));
		}
		ifs.close();
	}
	else
	{
		printf("\ninput error in filenames file");
	}
	return dataset;
}

void solve_dataset(int max_rep = 5)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(std::cin, filename);
	if (filename.size() == 0)
	{
		filename = "all.txt";
	}
	ofstream output;
	output.open("output.txt", ios::out);
	output << "solution methods for the CTOP \n";
	output << "filename,bestscore,score,cpu,gap\n";
	output.close();
	vector<Instance> set = read_dataset(filename);
	vector<Instance>::iterator it;
	for (it = set.begin(); it != set.end(); ++it)
	{
		it->result.resize(max_rep);
	}

	for (int rep = 0; rep < max_rep; ++rep)
	{
		//solve the dataset
		double avggap = 0.0;
		double avgscore = 0.0;
		for (it = set.begin(); it != set.end(); ++it)
		{
			Ins::MCTDTOPTW textfile = { it->path,it->filename };
			Ins instance(textfile);
			instance.read_time_independent_traveltime();
			instance.read_time_dependent_traveltime();
			instance.create_neighbourhood(textfile.path, textfile.name);
			//Aco acs(instance, 1, 3, 0.1, 20, 10000, 0.25, 0.05);
			//it->result[rep]=acs.solve(it->bestscore);
			Tabu tabu(instance, 10000,2);
			it->result[rep] = tabu.solve(it->bestknown);
			//cout << it->result[rep].sol << endl;
			//Ils ils(instance, 10000, 100, 20, 30);
			//it->result[rep] = ils.solve(it->bestscore);
			avggap += it->result[rep].gap;
			avgscore += it->result[rep].sol.score;
			cout << "name: " << it->filename << " best score: " << it->bestknown << " score: " << it->result[rep].sol.score << " cpu time: " << it->result[rep].time << " gap: " << it->result[rep].gap << endl;
			output.open("output.txt", ios::out | ios::app);
			output << it->filename << ";" << it->bestknown << ";" << it->result[rep].sol.score << ";" << it->result[rep].time << ";" << it->result[rep].gap << '\n';
			output.close();
		}
	}
	//calculate results over all replicates
	double avg_avg_gap = 0.0;
	double sd_avg_gap = 0.0;
	output.open("output.txt", ios::out | ios::app);
	output << "filename;avgscore;avgtime;avggap;sdgap\n";
	for (it = set.begin(); it != set.end(); ++it)
	{
		it->calculate_statistics();
		output<<it->filename<<";" << it->avgscore << ";" << it->avgtime<< ";"<<it->avggap<<";"<<it->sdgap << '\n';
		avg_avg_gap += it->avggap;
		
	}
	avg_avg_gap /= set.size();
	for (it = set.begin(); it != set.end(); ++it)
	{
		sd_avg_gap += pow(it->avggap-avg_avg_gap, 2);
	}
	sd_avg_gap = sqrt(sd_avg_gap / (set.size() - 1));
	output << "avg avg gap: " << avg_avg_gap<< " sd avg gap: "<< sd_avg_gap << '\n';
	output.close();
	cout << "avg avg gap is: " << avg_avg_gap<<" sd avg gap: "<<sd_avg_gap << endl;
}

void case_study(int max_rep = 5) 
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(std::cin, filename);
	if (filename.size() == 0)
	{
		filename = "case.txt";
	}
	ofstream output;
	output.open("output.txt", ios::out);
	output << "solution methods for the CTOP \n";
	output << "filename,bestscore,score,cpu,gap,removed\n";
	output.close();
	vector<Instance> dataset = read_dataset(filename);
	vector<Instance>::iterator it;
	for (it = dataset.begin(); it != dataset.end(); ++it)
	{
		it->result.resize(max_rep);
	}

	for (int rep = 0; rep < max_rep; ++rep)
	{
		//solve the dataset
		double avggap = 0.0;
		double avgscore = 0.0;
		for (it = dataset.begin(); it != dataset.end(); ++it)
		{
			Ins::MCTDTOPTW textfile = { it->path,it->filename };
			Ins instance(textfile);
			instance.read_time_independent_traveltime();
			instance.read_time_dependent_traveltime();
			instance.create_neighbourhood(textfile.path, textfile.name);
			//instance.alter_instance();
			Tabu tabu(instance, 10000, 2);
			it->result[rep] = tabu.solve(it->bestknown);
			//instance.unalter_instance();
			//cout << "after repair" << endl;
			//it->result[rep].removed=it->result[rep].sol.repair();
			cout << it->result[rep].sol << endl;
			avggap += it->result[rep].gap;
			avgscore += it->result[rep].sol.score;
			cout << "name: " << it->filename << " best score: " << it->bestknown << " score: " << it->result[rep].sol.score << " cpu time: " << it->result[rep].time << " gap: " << it->result[rep].gap<<" removed: "<< it->result[rep].removed << endl;
			output.open("output.txt", ios::out | ios::app);
			output << it->filename << ";" << it->bestknown << ";" << it->result[rep].sol.score << ";" << it->result[rep].time << ";" << it->result[rep].gap<< ";" << it->result[rep].removed << "\n";
			output.close();
		}
	}
	//calculate results over all replicates
	double globalgap = 0.0;
	output.open("output.txt", ios::out | ios::app);
	for (it = dataset.begin(); it != dataset.end(); ++it)
	{
		double avgscore = 0.0;
		double avgcpu = 0.0;
		double avgremoved = 0.0;
		for (int rep = 0; rep < max_rep; ++rep)
		{
			avgscore += it->result[rep].sol.score;
			avgcpu += it->result[rep].time;
			avgremoved += it->result[rep].removed;
		}
		avgscore /= max_rep;
		avgcpu /= max_rep;
		avgremoved /= max_rep;
		output << it->filename << ";" << avgscore << ";" << avgcpu << " ; " << avgremoved<< "\n";
		double avggap = (double(it->bestknown - avgscore) / it->bestknown) * 100;
		globalgap += avggap;
	}
	output.close();
	globalgap /= dataset.size();
	cout << "global avg gap is: " << globalgap << endl;

}

void debug_instance()
{
	Res resdebug;
	Ins::MCTDTOPTW textfile = { "..\\..\\datasets\\MCTDTOPTW\\" ,"100.3.1.2.txt" };
	Ins instance(textfile);
	instance.read_time_independent_traveltime();
	instance.read_time_dependent_traveltime();
	instance.create_neighbourhood(textfile.path, textfile.name);
	//Aco acs(instance, 1, 3, 0.01, 20, 10000, 0.25, 0.05);
	//resdataset.push_back(acs.solve());
	Tabu tabu(instance,10000,2);
	resdebug=tabu.solve(894);
	//Ils ils(instance, 10000, 100, 20, 30);
	//resdataset.push_back(ils.solve());
	cout << " best score: " << 894 << " score: "<< resdebug.sol.score << " cpu time: " << resdebug.time << " gap: " << resdebug.gap << endl;
}

void doe(int max_rep = 10)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(std::cin, filename);
	if (filename.size() == 0)
	{
		filename = "all.txt";
	}
	ofstream output;
	output.open("output.txt", ios::out);
	output << "DOE for TS \n";
	output << "Nnimax,umax,avg_gap,std_gap\n";
	output.close();
	vector<int>umax{2,4,6};
	vector<int>nimax{5000,10000,20000};
	vector<Instance> dataset = read_dataset(filename);
	vector<Instance>::iterator it;
	for (int par1 = 0; par1 < 3; ++par1)
	{
		for (int par2 = 0; par2 < 3; ++par2)
		{
			for (it = dataset.begin(); it != dataset.end(); ++it)
			{
				it->result.resize(max_rep);
			}
			for (int rep = 0; rep < max_rep; ++rep)
			{
				//solve the dataset
				double avggap = 0.0;
				double avgscore = 0.0;
				for (it = dataset.begin(); it != dataset.end(); ++it)
				{
					Ins::MCTDTOPTW textfile = { it->path,it->filename };
					Ins instance(textfile);
					instance.read_time_independent_traveltime();
					instance.read_time_dependent_traveltime();
					instance.create_neighbourhood(textfile.path, textfile.name);
					Tabu tabu(instance,nimax[par1], umax[par2]);
					it->result[rep] = tabu.solve(it->bestknown);
					avggap += it->result[rep].gap;
					avgscore += it->result[rep].sol.score;
					cout << "name: " << it->filename << " best score: " << it->bestknown << " score: " << it->result[rep].sol.score << " cpu time: " << it->result[rep].time << " gap: " << it->result[rep].gap << endl;
				}//end it
			}//end rep
			//calculate results over all replicates
			double globalgap = 0.0;
			double globalgap_sq = 0.0;
			for (it = dataset.begin(); it != dataset.end(); ++it)
			{
				double avgscore = 0.0;
				double avgcpu = 0.0;
				vector<double> gaps(max_rep);

				for (int rep = 0; rep < max_rep; ++rep)
				{
					avgscore += it->result[rep].sol.score;
					avgcpu += it->result[rep].time;
					gaps[rep] = it->result[rep].gap;
				}

				avgscore /= max_rep;
				avgcpu /= max_rep;

				// Compute avg gap per instance
				double avggap = (double(it->bestknown - avgscore) / it->bestknown) * 100.0;
				globalgap += avggap;

				// Compute standard deviation of gap
				double gap_sum = accumulate(gaps.begin(), gaps.end(), 0.0);
				double gap_mean = gap_sum / max_rep;

				double sq_sum = 0.0;
				for (double g : gaps)
					sq_sum += (g - gap_mean) * (g - gap_mean);

				double gap_stdev = sqrt(sq_sum / max_rep);
				globalgap_sq += gap_stdev;
				output.open("output.txt", ios::out | ios::app);
				output << it->filename << ","<< nimax[par1] << "," << umax[par2] << ","<< gap_mean <<","<< gap_stdev << "\n";

			}
			globalgap /= dataset.size();
			globalgap_sq /= dataset.size();

			cout<<"nimax: " << nimax[par1]<<"umax: " << umax[par2] << " avg gap is: " << globalgap <<"sd gap is: "<<globalgap_sq << endl;
			output.open("output.txt", ios::out | ios::app);
			output << nimax[par1] << ";" << umax[par2] << ";" << globalgap << ";" << globalgap_sq << "\n";
			output.close();
		}//end par 1
	}//end par2
}//end doe

void debug_ctop()
{
	Res res;
	//Ins::CTOP textfile = {"..\\..\\datasets\\CTOP\\LargeScale CTOP\\set2\\","b89.txt"};
	Ins::CTOP textfile = {"..\\..\\datasets\\CTOP\\DatasetsCTOP\\2set\\","b5.txt"};
	Ins instance(textfile);
	instance.create_neighbourhood(textfile.path, textfile.name);
	//Ils ils(instance, 10000, 100, 20, 30);
	//res = ils.solve();
	Tabu tabu(instance, 10000,2);
	res = tabu.solve(139);
	//Aco acs(instance, 1, 2, 0.01, 20, 10000, 0.25, 0.05);
	//res=acs.solve();
	cout << res.sol.score << " cpu time: " << res.time << endl;
}

void ctop_gap(int max_rep=5)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(std::cin, filename);
	if (filename.size() == 0)
	{
		filename = "set5.txt";
	}
	ofstream output;
	output.open("output.txt", ios::out);
	output << "solution methods for the CTOP: " << filename << "\n" << endl;
	output << "filename,bestscore,score,cpu,gap\n";
	output.close();
	vector<Instance> dataset = read_dataset(filename);
	vector<Instance>::iterator it;
	for (it = dataset.begin(); it != dataset.end(); ++it)
	{
		it->result.resize(max_rep);
	}
	
	for (int rep = 0; rep < max_rep; ++rep)
	{
		//solve the dataset
		double avggap = 0.0;
		double avgscore = 0.0;
		for (it = dataset.begin(); it != dataset.end(); ++it)
		{
			Ins::CTOP textfile = { it->path,it->filename };
			Ins instance(textfile);
			instance.create_neighbourhood(textfile.path, textfile.name);
			//Aco acs(instance, 1,1, 0.01, 20, 10000, 0.25, 0.05);
			//it->result[rep] = acs.solve(it->bestscore);
			//Ils ils(instance, 10000,100,2,3);
			//it->result[rep] = ils.solve(it->bestscore);
			Tabu tabu(instance, 10000,2);
			it->result[rep] = tabu.solve(it->bestknown);
			avggap += it->result[rep].gap;
			avgscore += it->result[rep].sol.score;
			cout << "name: " << it->filename << " best score: " << it->bestknown<<" " << tabu.name << " score: " << it->result[rep].sol.score << " cpu time: " << it->result[rep].time << " gap: " << it->result[rep].gap << endl;
			output.open("output.txt", ios::out | ios::app);
			output << it->filename << ";" << it->bestknown << ";" << it->result[rep].sol.score << ";" << it->result[rep].time << ";" << it->result[rep].gap << "\n";
			output.close();
		}
		//calculate dataset performance
		avggap /= dataset.size();
		avgscore /= dataset.size();
		cout << "average gap of " << filename << " is: " << avggap << " avg score: " << avgscore << endl;
	}
	//calculate results over all replicates
	double globalgap = 0.0;
	for (it = dataset.begin(); it != dataset.end(); ++it)
	{
		double avgscore = 0.0;
		for (int rep = 0; rep < max_rep; ++rep)
		{
			avgscore += it->result[rep].sol.score;
		}
		avgscore /= max_rep;
		double avggap= (double(it->bestknown - avgscore) / it->bestknown) * 100;
		globalgap += avggap;
	}
	globalgap /= dataset.size();
	cout << "global avg gap is: " << globalgap << endl;
}//end ctop_gap

int main()
{
	//create_case_dataset();
	//Graph bemobile(425479, 519915);
	//debug_instance();
	//debug_ctop(5);
	solve_dataset(1);
	//ctop_gap(5);
	//doe(10);
	//case_study(10);


}
