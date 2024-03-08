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
