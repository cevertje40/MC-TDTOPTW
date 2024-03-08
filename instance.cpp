#include "instance.h"

Ins::Ins(string name)
{
	char filepath[125] = "..\\..\\datasets\\MCTDTOPTW\\";
	strcat(filepath, name);
	FILE* file;
	fopen_s(&file, filepath, "r");
	if (!file)
	{
		printf("\ninput error in instance file");
	}
	fscanf(file, "%d", &maxvertices);
	fscanf(file, "%d", &paths);
	score.resize(maxvertices, 0);
	wei.resize(maxvertices, 0);
	vol.resize(maxvertices, 0);
	LTW.resize(maxvertices, vector<double>(paths, 0));
	UTW.resize(maxvertices, vector<double>(paths, 0));
	serv.resize(maxvertices);
	T_max.resize(paths);
	wei_max.resize(paths);
	vol_max.resize(paths);
	fscanf(file, "%lf %lf %lf", &brk_dur, &brk_start, &brk_end);
	for (int p = 0; p < paths; ++p)
	{
		fscanf(file, "%lf %lf %lf", &T_max[p], &wei_max[p], &vol_max[p]);
	}
	for (int i = 0; i < maxvertices; ++i)
	{
		int temp;
		fscanf_s(file, "%d %d %lf %lf %lf", &temp, &score[i], &serv[i], &wei[i], &vol[i]);
		for (int p = 0; p < paths; ++p)
		{
			fscanf_s(file, "%lf %lf", &LTW[i][p], &UTW[i][p]);
		}
	}
	fclose(file);
}
