// MC-TDTOPTW.cpp


#include "solver.h"
#include "graph.h"

using namespace std;
int main()
{
    Graph bemobile(425479, 519915, 56);
    //cout << bemobile.dijkstra_independent(41, 18467) << endl;
    //cout << bemobile.dijkstra_dependent(41,18467,6) << endl;
    string name = "20.txt";
    Ins instance(name);
    //instance.construct_time_independent_traveltime(bemobile);
    instance.construct_time_dependent_traveltime(bemobile);
    Solver solver;
    solver.aco(instance,1,3,0.01,20,10000);
}
