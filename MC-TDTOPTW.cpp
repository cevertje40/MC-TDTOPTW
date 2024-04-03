// MC-TDTOPTW.cpp

#include "solver.h"

using namespace std;


int main()
{
    //Graph bemobile(425479, 519915);
    //cout << bemobile.dijkstra_independent(41, 18467) << endl;
    //cout << bemobile.dijkstra_dependent(41,18467,6) << endl;

    //create datasets
    //create or read neighbourhood
    //instance.create_neighbourhood(45);
    string name = "20.txt";
    Ins instance(name);
    instance.read_neighbourhood();

    //read instance info and travel time
   
    //aco,ils class
    //construction aco
    //local search moves
    //test class
    
    Aco acs(instance, 1,3,0.01,20,10000);
    acs.solve();
}
