// MC-TDTOPTW.cpp

#include "solver.h"

using namespace std;
int main()
{
    //Graph bemobile(425479, 519915);
    //cout << bemobile.dijkstra_independent(41, 18467) << endl;
    //cout << bemobile.dijkstra_dependent(41,18467,6) << endl;
    string name = "20.txt";
    Ins instance(name);
   
    Solver solver;
    solver.aco(instance,1,3,0.01,20,10000);
}
