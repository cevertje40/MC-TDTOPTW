// MC-TDTOPTW.cpp

#include "solver.h"

using namespace std;
int main()
{
    //Graph bemobile(425479, 519915);
    //cout << bemobile.dijkstra_independent(41, 18467) << endl;
    //cout << bemobile.dijkstra_dependent(41,18467,6) << endl;

    //read instance info and travel time
    string name = "20.txt";
    Ins instance(name);

    //create or read neighbourhood
    instance.create_neighbourhood(45);
    instance.read_neighbourhood();
   
    //read neigbhorhood
    //create datasets
    //aco,ils class
    //construction aco
    //local search moves
    //test class

    Solver solver;
    solver.aco(instance,1,3,0.01,20,10000);
}
