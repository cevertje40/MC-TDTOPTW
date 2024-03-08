// MC-TDTOPTW.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "solver.h"

using namespace std;
int main()
{
    Ins instance;
    string name = "8.txt";
    instance.read_instance(name);
    Solver solver;
    solver.aco(instance,1,3,0.01,20,10000);
    cout << "Hello World" << endl;;
}
