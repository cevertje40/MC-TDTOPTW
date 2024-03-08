// MC-TDTOPTW.cpp


#include "solver.h"
#include "graph.h"

using namespace std;
int main()
{
    string name = "8.txt";
    Ins instance(name);
    Solver solver;
    solver.aco(instance,1,3,0.01,20,10000);
    cout << "Hello World" << endl;;
}
