// MC-TDTOPTW.cpp

#include "solver.h"
using namespace std;

class Instance
{
public:
	string path;
	string filename;
	vector<Res> result;
	int bestknown{ 0 };
	double avgscore{ 0.0 };
	double sdscore{ 0.0 };
	double avggap{ 0.0 };
	double sdgap{ 0.0 };
	double avgtime{ 0.0 };
	Instance(const string& path, const string& filename, int bestknown) : path(path), filename(filename), bestknown(bestknown) {}
	void calculate_statistics()
	{
		avgscore = 0.0;
		sdscore = 0.0;
		avggap = 0.0;
		sdgap = 0.0;
		avgtime = 0.0;
		const int n = (int)result.size();

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
		sdscore = sqrt(sdscore / (n - 1));
		sdgap = sqrt(sdgap / (n - 1));
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
			//instance.create_neighbourhood_simple(textfile.path, textfile.name);
			instance.create_neighbourhood(textfile.path, textfile.name);
			//Aco acs(instance, 1, 3, 0.1, 20, 10000, 0.25, 0.05);
			//it->result[rep]=acs.solve(it->bestscore);
			

			double maxtime=9999;
			/*
			if (instance.maxvertices == 20)
				maxtime = 1;
			else if (instance.maxvertices == 50)
				maxtime = 2.0;
			else
				maxtime = 5.0;
			*/

			Tabu tabu(instance, 20000, 2, 0.1, 0.9, 10);
			it->result[rep] = tabu.solve(it->bestknown,maxtime,-1,7);
			/*
			forced_crit = -1 for adaptive
				forced_crit = 0 for SCORE
				forced_crit = 1 for TIME
				forced_crit = 2 for VOLUME
				forced_crit = 3 for WEIGHT
			*/
			//it->result[rep].sol.write_to_cplex(textfile.name);
			//cout << it->result[rep].sol << endl;
			//Ils ils(instance, 10000, 100, 20, 30);
			//it->result[rep] = ils.solve(it->bestscore);
			//Alns alns(instance,5000,0.9997,150,500,8.0,4.0,1.0,100.0,0.0001,0.15,0.8);
			//it->result[rep] = alns.solve(it->bestknown,maxtime);

			avggap += it->result[rep].gap;
			avgscore += it->result[rep].sol.score;
			cout << "name: " << it->filename << " best score: " << it->bestknown << " score: " << it->result[rep].sol.score << " cpu time: " << it->result[rep].time << " gap: " << it->result[rep].gap << endl;
			output.open("output.txt", ios::out | ios::app);
			output << it->filename << ";" << it->bestknown << ";" << it->result[rep].sol.score << ";" << it->result[rep].time << ";" << it->result[rep].gap << '\n';
			output.close();
		}
	}
	//calculate results over all replicates
	double avg_cpu = 0.0;
	double avg_avg_gap = 0.0;
	double sd_avg_gap = 0.0;
	output.open("output.txt", ios::out | ios::app);
	output << "filename;avgscore;avgtime;avggap;sdgap\n";
	for (it = set.begin(); it != set.end(); ++it)
	{
		it->calculate_statistics();
		output<<it->filename<<";" << it->avgscore << ";" << it->avgtime<< ";"<<it->avggap<<";"<<it->sdgap << '\n';
		avg_avg_gap += it->avggap;
		avg_cpu += it->avgtime;
		
	}
	avg_avg_gap /= set.size();
	avg_cpu /= set.size();
	for (it = set.begin(); it != set.end(); ++it)
	{
		sd_avg_gap += pow(it->avggap-avg_avg_gap, 2);
	}
	sd_avg_gap = sqrt(sd_avg_gap / (set.size() - 1));
	output << "avg avg gap: " << avg_avg_gap<< " sd avg gap: "<< sd_avg_gap<< " avg cpu: " << avg_cpu << '\n';
	output.close();
	cout << "avg avg gap is: " << avg_avg_gap<<" sd avg gap: "<<sd_avg_gap << " avg cpu: " << avg_cpu << endl;
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
	output << "case study results \n";
	output << "filename;bestscore;score;cpu;gap;total_served_customers;total_removed_customers;avg_score_per_served_customer;avg_route_duration;avg_route_time_utilization;avg_weight_utilization;avg_volume_utilization;avg_waiting_time_per_customer;avg_break_start_time;avg_break_position_norm;pct_break_at_end_depot;avg_removed_score;avg_removed_tw_width;avg_removed_service;avg_removed_weight;avg_removed_volume;avg_removed_depot_tt;avg_removed_position\n";
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
		for (it = dataset.begin(); it != dataset.end(); ++it)
		{
			Ins::MCTDTOPTW textfile = { it->path,it->filename };
			Ins instance(textfile);
			instance.read_time_independent_traveltime();
			instance.read_time_dependent_traveltime();
			instance.create_neighbourhood(textfile.path, textfile.name);
			//instance.alter_instance();
			Tabu tabu(instance, 10000, 2, 0.1, 0.9, 10);
			it->result[rep] = tabu.solve(it->bestknown);
			//Alns alns(instance,5000,0.9997,150,200,8.0,4.0,1.0,100.0,0.0001,0.15,0.8);
			//it->result[rep] = alns.solve(it->bestknown);
			//instance.unalter_instance();
			//cout << "after repair" << endl;
		
			it->result[rep].stats = it->result[rep].sol.repair_and_collect_stats();

			cout << "name: " << it->filename
				<< " best score: " << it->bestknown
				<< " score: " << it->result[rep].sol.score
				<< " cpu time: " << it->result[rep].time
				<< " gap: " << it->result[rep].gap
				<< " served: " << it->result[rep].stats.total_served_customers
				<< " removed: " << it->result[rep].stats.total_removed_customers
				<< " avg score/ served customer: " << it->result[rep].stats.avg_score_per_served_customer
				<< " avg route duration: " << it->result[rep].stats.avg_route_duration
				<< " avg route time utilization: " << it->result[rep].stats.avg_route_time_utilization
				<< " avg weight utilization: " << it->result[rep].stats.avg_weight_utilization
				<< " avg volume utilization: " << it->result[rep].stats.avg_volume_utilization
				<< " avg waiting time/customer: " << it->result[rep].stats.avg_waiting_time_per_customer
				<< " avg break start time: " << it->result[rep].stats.avg_break_start_time
				<< " avg break position norm: " << it->result[rep].stats.avg_break_position_norm
				<< " pct break at end depot: " << it->result[rep].stats.pct_break_at_end_depot
				<< " avg removed score: " << it->result[rep].stats.avg_removed_score
				<< " avg removed tw width: " << it->result[rep].stats.avg_removed_tw_width
				<< " avg removed service: " << it->result[rep].stats.avg_removed_service
				<< " avg removed weight: " << it->result[rep].stats.avg_removed_weight
				<< " avg removed volume: " << it->result[rep].stats.avg_removed_volume
				<< " avg removed depot tt: " << it->result[rep].stats.avg_removed_depot_tt
				<< " avg removed position: " << it->result[rep].stats.avg_removed_position
				<< endl;
			
			cout << it->result[rep].sol << endl;
			output.open("output.txt", ios::out | ios::app);
			output << it->filename << ";"
				<< it->bestknown << ";"
				<< it->result[rep].sol.score << ";"
				<< it->result[rep].time << ";"
				<< it->result[rep].gap << ";"
				<< it->result[rep].stats.total_served_customers << ";"
				<< it->result[rep].stats.total_removed_customers << ";"
				<< it->result[rep].stats.avg_score_per_served_customer << ";"
				<< it->result[rep].stats.avg_route_duration << ";"
				<< it->result[rep].stats.avg_route_time_utilization << ";"
				<< it->result[rep].stats.avg_weight_utilization << ";"
				<< it->result[rep].stats.avg_volume_utilization << ";"
				<< it->result[rep].stats.avg_waiting_time_per_customer << ";"
				<< it->result[rep].stats.avg_break_start_time << ";"
				<< it->result[rep].stats.avg_break_position_norm << ";"
				<< it->result[rep].stats.pct_break_at_end_depot << ";"
				<< it->result[rep].stats.avg_removed_score << ";"
				<< it->result[rep].stats.avg_removed_tw_width << ";"
				<< it->result[rep].stats.avg_removed_service << ";"
				<< it->result[rep].stats.avg_removed_weight << ";"
				<< it->result[rep].stats.avg_removed_volume << ";"
				<< it->result[rep].stats.avg_removed_depot_tt << ";"
				<< it->result[rep].stats.avg_removed_position << ";"
				<< "\n";
			output.close();
		}
	}
	//calculate results over all replicates
	double globalgap = 0.0;
	output.open("output.txt", ios::out | ios::app);
	output << "\nsummary over replications\n";
	output << "filename;avg_score;avg_cpu;avg_removed;avg_served;avg_score_per_customer;avg_route_duration;avg_route_util;avg_weight_util;avg_volume_util;avg_waiting_per_customer;avg_break_start;avg_break_pos;avg_pct_break_end_depot;avg_removed_score;avg_removed_tw;avg_removed_service;avg_removed_weight;avg_removed_volume;avg_removed_depot_tt;avg_removed_position\n";

	for (it = dataset.begin(); it != dataset.end(); ++it)
	{
		double avgscore = 0.0;
		double avgcpu = 0.0;
		double avgremoved = 0.0;
		double avgserved = 0.0;
		double avgscorepercust = 0.0;
		double avgroute_duration = 0.0;
		double avgrouteutil = 0.0;
		double avgweightutil = 0.0;
		double avgvolumeutil = 0.0;
		double avgwait = 0.0;
		double avgbreakstart = 0.0;
		double avgbreakpos = 0.0;
		double avgpctbreakend = 0.0;
		double avgremovedscore = 0.0;
		double avgremovedtw = 0.0;
		double avgremovedservice = 0.0;
		double avgremovedweight = 0.0;
		double avgremovedvolume = 0.0;
		double avgremoveddepot = 0.0;
		double avgremovedposition = 0.0;

		for (int rep = 0; rep < max_rep; ++rep)
		{
			avgscore += it->result[rep].sol.score;
			avgcpu += it->result[rep].time;
			avgremoved += it->result[rep].stats.total_removed_customers;
			avgserved += it->result[rep].stats.total_served_customers;
			avgscorepercust += it->result[rep].stats.avg_score_per_served_customer;
			avgroute_duration += it->result[rep].stats.avg_route_duration;
			avgrouteutil += it->result[rep].stats.avg_route_time_utilization;
			avgweightutil += it->result[rep].stats.avg_weight_utilization;
			avgvolumeutil += it->result[rep].stats.avg_volume_utilization;
			avgwait += it->result[rep].stats.avg_waiting_time_per_customer;
			avgbreakstart += it->result[rep].stats.avg_break_start_time;
			avgbreakpos += it->result[rep].stats.avg_break_position_norm;
			avgpctbreakend += it->result[rep].stats.pct_break_at_end_depot;
			avgremovedscore += it->result[rep].stats.avg_removed_score;
			avgremovedtw += it->result[rep].stats.avg_removed_tw_width;
			avgremovedservice += it->result[rep].stats.avg_removed_service;
			avgremovedweight += it->result[rep].stats.avg_removed_weight;
			avgremovedvolume += it->result[rep].stats.avg_removed_volume;
			avgremoveddepot += it->result[rep].stats.avg_removed_depot_tt;
			avgremovedposition += it->result[rep].stats.avg_removed_position;
		}

		avgscore /= max_rep;
		avgcpu /= max_rep;
		avgremoved /= max_rep;
		avgserved /= max_rep;
		avgscorepercust /= max_rep;
		avgroute_duration /= max_rep;
		avgrouteutil /= max_rep;
		avgweightutil /= max_rep;
		avgvolumeutil /= max_rep;
		avgwait /= max_rep;
		avgbreakstart /= max_rep;
		avgbreakpos /= max_rep;
		avgpctbreakend /= max_rep;
		avgremovedscore /= max_rep;
		avgremovedtw /= max_rep;
		avgremovedservice /= max_rep;
		avgremovedweight /= max_rep;
		avgremovedvolume /= max_rep;
		avgremoveddepot /= max_rep;
		avgremovedposition /= max_rep;

		output << it->filename << ";"
			<< avgscore << ";"
			<< avgcpu << ";"
			<< avgremoved << ";"
			<< avgserved << ";"
			<< avgscorepercust << ";"
			<< avgroute_duration << ";"
			<< avgrouteutil << ";"
			<< avgweightutil << ";"
			<< avgvolumeutil << ";"
			<< avgwait << ";"
			<< avgbreakstart << ";"
			<< avgbreakpos << ";"
			<< avgpctbreakend << ";"
			<< avgremovedscore << ";"
			<< avgremovedtw << ";"
			<< avgremovedservice << ";"
			<< avgremovedweight << ";"
			<< avgremovedvolume << ";"
			<< avgremoveddepot << ";"
			<< avgremovedposition
			<< "\n";

		double avggap = (double(it->bestknown - avgscore) / it->bestknown) * 100.0;
		globalgap += avggap;
	}

	output.close();
	globalgap /= dataset.size();
	cout << "global avg gap is: " << globalgap << endl;
}

void debug_instance()
{
	Res resdebug;
	Ins::MCTDTOPTW textfile = { "..\\..\\datasets\\MCTDTOPTW\\" ,"20.1.3.2.txt" };
	Ins instance(textfile);
	instance.read_time_independent_traveltime();
	instance.read_time_dependent_traveltime();
	instance.create_neighbourhood(textfile.path, textfile.name);
	//Aco acs(instance, 1, 3, 0.01, 20, 10000, 0.25, 0.05);
	//resdataset.push_back(acs.solve());
	Tabu tabu(instance, 10000, 2, 0.1, 0.8, 5);
	resdebug=tabu.solve(416);
	resdebug.sol.write_to_cplex(textfile.name);
	//Ils ils(instance, 10000, 100, 20, 30);
	//resdataset.push_back(ils.solve());
	//Alns alns(instance, 5000, 0.9997, 150, 1500, 8.0, 4.0, 1.0, 100.0, 0.0001, 0.15, 0.8);
	//resdebug = alns.solve(250);
	cout << " best score: " << 416 << " score: "<< resdebug.sol.score << " cpu time: " << resdebug.time << " gap: " << resdebug.gap << endl;
}

void doe(int max_rep = 10)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	std::string filename;
	std::getline(cin, filename);
	if (filename.empty()) filename = "all.txt";

	// Parameter grids
	std::vector<int> umax{ 2, 4, 6 };
	std::vector<int> nimax{ 5000, 10000, 20000 };

	// Read dataset
	std::vector<Instance> dataset = read_dataset(filename);

	// Output files (separate, unambiguous)
	{
		std::ofstream instout("doe1_instance_means.csv", std::ios::out);
		instout << "Nnimax;umax;instance;mean_gap;sd_gap;mean_time\n";
	}
	{
		std::ofstream overout("doe1_overall.csv", std::ios::out);
		overout << "Nnimax;umax;mean_gap_over_instances;sd_gap_over_instances;mean_time_over_instances\n";
	}

	// Helpers
	auto mean_of = [](const vector<double>& v) 
	{
		return v.empty() ? 0.0 : accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
	};
	auto sd_sample_of = [&](const vector<double>& v) 
	{
		if (v.size() < 2) return 0.0;
		double m = mean_of(v), acc = 0.0;
		for (double x : v) acc += (x - m) * (x - m);
		return sqrt(acc / static_cast<double>(v.size() - 1));
	};

	for (int p1 = 0; p1 < static_cast<int>(nimax.size()); ++p1) 
	{
		for (int p2 = 0; p2 < static_cast<int>(umax.size()); ++p2) 
		{
			// Per-instance aggregates (means over reps)
			std::vector<double> per_inst_mean_gap;  per_inst_mean_gap.reserve(dataset.size());
			std::vector<double> per_inst_mean_time; per_inst_mean_time.reserve(dataset.size());

			// Prepare storage for results
			for (auto& inst : dataset) inst.result.resize(max_rep);

			// Run all reps for all instances for this (nimax, umax)
			for (int rep = 0; rep < max_rep; ++rep) 
			{
				for (auto& inst : dataset) 
				{
					Ins::MCTDTOPTW textfile{ inst.path, inst.filename };
					Ins instance(textfile);
					instance.read_time_independent_traveltime();
					instance.read_time_dependent_traveltime();
					instance.create_neighbourhood(textfile.path, textfile.name);

					// Keep the epsilon, theta, gamma fixed as in your DOE1
					Tabu tabu(instance, nimax[p1], umax[p2], 0.1, 0.8, 5);
					inst.result[rep] = tabu.solve(inst.bestknown);

					cout << "name: " << inst.filename<< " best score: " << inst.bestknown<< " score: " << inst.result[rep].sol.score<< " cpu time: " << inst.result[rep].time<< " gap: " << inst.result[rep].gap << endl;
				}
			}

			// Write per-instance rows + accumulate for overall
			{
				ofstream instout("doe1_instance_means.csv", std::ios::app);
				for (auto& inst : dataset) 
				{
					std::vector<double> gaps(max_rep), times(max_rep);
					for (int rep = 0; rep < max_rep; ++rep) 
					{
						gaps[rep] = inst.result[rep].gap;   // already (% gap)
						times[rep] = inst.result[rep].time;  // seconds
					}

					double mean_gap = mean_of(gaps);
					double sd_gap = sd_sample_of(gaps);     // sample SD across reps
					double mean_time = mean_of(times);

					per_inst_mean_gap.push_back(mean_gap);
					per_inst_mean_time.push_back(mean_time);

					instout << nimax[p1] << ';'<< umax[p2] << ';'<< inst.filename << ';'<< mean_gap << ';'<< sd_gap << ';'<< mean_time << "\n";
				}
			}

			// Overall across instances (of per-instance means)
			double mean_gap_over_inst = mean_of(per_inst_mean_gap);
			double sd_gap_over_inst = sd_sample_of(per_inst_mean_gap); // variability across instances
			double mean_time_over_inst = mean_of(per_inst_mean_time);
			{
				ofstream overout("doe1_overall.csv", std::ios::app);
				overout << nimax[p1] << ';'	<< umax[p2] << ';'<< mean_gap_over_inst << ';'<< sd_gap_over_inst << ';'<< mean_time_over_inst << "\n";
			}

			std::cout << "Nnimax: " << nimax[p1]<< " umax: " << umax[p2]<< " | mean_gap: " << mean_gap_over_inst<< " sd_gap: " << sd_gap_over_inst<< " mean_time: " << mean_time_over_inst<< std::endl;
		}
	}
}

void doe2(int max_rep = 5)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(std::cin, filename);
	if (filename.size() == 0)
	{
		filename = "all.txt";
	}
	vector<double>epsilon{0.05,0.10,0.2};
	vector<double>theta{0.7,0.8,0.9};
	vector<double>gamma{2.5,5.0,10.0};
	vector<Instance> dataset = read_dataset(filename);
	std::ofstream instout("doe2_instance_means.csv", std::ios::out);
	instout << "epsilon;theta;gamma;instance;mean_gap;sd_gap;mean_time\n";
	instout.close();
	std::ofstream overout("doe2_overall.csv", std::ios::out);
	overout << "epsilon;theta;gamma;mean_gap_over_instances;sd_gap_over_instances;mean_time_over_instances\n";
	overout.close();
	for (int par1 = 0; par1 < 3; ++par1)
	{
		for (int par2 = 0; par2 < 3; ++par2)
		{
			for (int par3 = 0; par3 < 3; ++par3)
			{
				// ---- Aggregate over reps: per-instance means (then overall across instances)
				std::vector<double> per_inst_mean_gap;
				std::vector<double> per_inst_mean_time;
				per_inst_mean_gap.reserve(dataset.size());
				per_inst_mean_time.reserve(dataset.size());
				std::ofstream instout("doe2_instance_means.csv", std::ios::app);
				for (auto& inst : dataset) 
				{
					std::vector<double> gaps(max_rep), times(max_rep);
					inst.result.resize(max_rep);
					for (int rep = 0; rep < max_rep; ++rep) 
					{
						Ins::MCTDTOPTW textfile = {inst.path,inst.filename};
						Ins instance(textfile); instance.read_time_independent_traveltime();
						instance.read_time_dependent_traveltime();
						instance.create_neighbourhood(textfile.path, textfile.name);
						Tabu tabu(instance, 10000, 2, epsilon[par1], theta[par2], gamma[par3]);
						inst.result[rep] = tabu.solve(inst.bestknown);
						gaps[rep] = inst.result[rep].gap;
						times[rep] = inst.result[rep].time;
						cout << "name: " << inst.filename << " best score: " << inst.bestknown << " score: " << inst.result[rep].sol.score << " cpu time: " << inst.result[rep].time << " gap: " << inst.result[rep].gap << endl;
					}
					// per-instance mean gap and time
					double mean_gap = std::accumulate(gaps.begin(), gaps.end(), 0.0) / max_rep;
					double mean_time = std::accumulate(times.begin(), times.end(), 0.0) / max_rep;

					// per-instance sample stdev of gap across reps
					double var_gap = 0.0;
					if (max_rep > 1) 
					{
						for (double g : gaps) var_gap += (g - mean_gap) * (g - mean_gap);
						var_gap /= (max_rep - 1);
					}
					double sd_gap = std::sqrt(var_gap);

					per_inst_mean_gap.push_back(mean_gap);
					per_inst_mean_time.push_back(mean_time);

					// write per-instance row (now unambiguous)
					
					instout << epsilon[par1] << ';' << theta[par2] << ';' << gamma[par3] << ';'<< inst.filename << ';' << mean_gap << ';' << sd_gap << ';' << mean_time << "\n";
					
				}
				instout.close();
				// overall mean across instances (of per-instance means)
				auto mean_of = [](const std::vector<double>& v) 
				{
					return std::accumulate(v.begin(), v.end(), 0.0) / std::max<size_t>(1, v.size());
				};
				auto sd_of = [&](const std::vector<double>& v) 
				{
					double m = mean_of(v), acc = 0.0;
					if (v.size() > 1) 
					{
						for (double x : v) acc += (x - m) * (x - m);
						acc /= (v.size() - 1);  // sample stdev across instances
					}
					return std::sqrt(acc);
				};

				double mean_gap_over_inst = mean_of(per_inst_mean_gap);
				double sd_gap_over_inst = sd_of(per_inst_mean_gap);
				double mean_time_over_inst = mean_of(per_inst_mean_time);

				// write 1 overall row per config to a separate file
				std::ofstream overout("doe2_overall.csv", std::ios::app);
				overout << epsilon[par1] << ';' << theta[par2] << ';' << gamma[par3] << ';'	<< mean_gap_over_inst << ';' << sd_gap_over_inst << ';' << mean_time_over_inst << "\n";
				overout.close();

				std::cout << "epsilon: " << epsilon[par1]
					<< " theta: " << theta[par2]
					<< " gamma: " << gamma[par3]
					<< " | mean_gap: " << mean_gap_over_inst
					<< " sd_gap: " << sd_gap_over_inst
					<< " mean_time: " << mean_time_over_inst << std::endl;
			}//end par3
		}//end par2
	}//end par1
}//end doe 2

void debug_ctop()
{
	Res res;
	Ins::CTOP textfile = {"..\\..\\datasets\\CTOP\\LargeScale CTOP\\set2\\","b2.txt"};
	//Ins::CTOP textfile = {"..\\..\\datasets\\CTOP\\DatasetsCTOP\\2set\\","b5.txt"};
	Ins instance(textfile);
	instance.create_neighbourhood(textfile.path, textfile.name);
	//Ils ils(instance, 10000, 100, 20, 30);
	//res = ils.solve();
	Tabu tabu(instance, 10000, 2, 0.1, 0.8, 5);
	res = tabu.solve(139);
	res.sol.write_to_file(textfile.name);
	//Aco acs(instance, 1, 2, 0.01, 20, 10000, 0.25, 0.05);
	//res=acs.solve();
	cout << res.sol.score << " cpu time: " << res.time << endl;
}

void ctop_gap(int max_rep=5)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(cin, filename);
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
			Tabu tabu(instance, 10000, 2, 0.1, 0.8, 5);
			it->result[rep] = tabu.solve(it->bestknown);
			avggap += it->result[rep].gap;
			avgscore += it->result[rep].sol.score;
			cout << "name: " << it->filename << " best score: " << it->bestknown<<" " << tabu.name << " score: " << it->result[rep].sol.score << " cpu time: " << it->result[rep].time << " gap: " << it->result[rep].gap << endl;
			output.open("output.txt", ios::out | ios::app);
			output << it->filename << ";" << it->bestknown << ";" << it->result[rep].sol.score << ";" << it->result[rep].time << ";" << it->result[rep].gap << "\n";
			output.close();
			//it->result[rep].sol.write_to_file(it->filename);
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

void ctop_optimal(int max_rep = 1)
{
	cout << fixed << setprecision(2) << "enter name of dataset" << endl;
	string filename;
	getline(cin, filename);
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
			Ins::KnownOptimalCTOP textfile = { it->path,it->filename };
			Ins instance(textfile);
			instance.create_neighbourhood(textfile.path, textfile.name);
			//Aco acs(instance, 1,1, 0.01, 20, 10000, 0.25, 0.05);
			//it->result[rep] = acs.solve(it->bestscore);
			//Ils ils(instance, 10000,100,2,3);
			//it->result[rep] = ils.solve(it->bestscore);
			Tabu tabu(instance, 10000,2, 0.1, 0.8, 5);
			it->result[rep] = tabu.solve(it->bestknown);
			Sol planted_sol(instance);
			planted_sol.read_from_file(textfile.name);
			it->bestknown = planted_sol.score;
			it->result[rep].gap = (double(it->bestknown - it->result[rep].sol.score) / it->bestknown) * 100;
			avggap += it->result[rep].gap;
			avgscore += it->result[rep].sol.score;
			cout << "name: " << it->filename << " planted score: " << it->bestknown << " " << tabu.name << " score: " << it->result[rep].sol.score << " cpu time: " << it->result[rep].time << " gap: " << it->result[rep].gap << endl;
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
		double avggap = (double(it->bestknown - avgscore) / it->bestknown) * 100;
		globalgap += avggap;
	}
	globalgap /= dataset.size();
	cout << "global avg gap is: " << globalgap << endl;

}

int main()
{
	//create_case_dataset();
	//Graph bemobile(425479, 519915);
	//debug_instance();
	//debug_ctop();
	solve_dataset(5);
	//ctop_gap();
	//doe(5);
	//doe2(5);
	//case_study(10);
	//ctop_optimal(5);

}
