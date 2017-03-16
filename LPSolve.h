#ifndef LPSOLVE_H
#define LPSOLVE_H
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include "algebra.h"
#include "interval.h"
using namespace std;
interval solveInterval(int varIdx, Box &box, deque<FastLinear>&ineqs, deque<FastLinear>&eqs, int numVars);
#endif 