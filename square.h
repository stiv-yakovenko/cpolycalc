#ifndef CPOLYCALC_SQUARE_H
#define CPOLYCALC_SQUARE_H

#include "Reduced.h"
#include <vector>
struct Replacement {
	Replacement() {
	}
	int x, y, z;
	double k;
};
void fillSquaresAndRoots(Reduced&reduced);
bool findSqSubst(std::map<std::pair<int, int>, double> &squares, int y, Replacement &repl);
map<int, int> calcPowers(const multiset<int> &vars);
void replaceVar(Polynome &poly, int v, int y, const multiset<int> &key,
	multiset<int> replacement, double k);
multiset<int> genMonomePowerKey(int x, int n);
bool isPower(Polynome &poly, int n, Replacement &rep);
int calcNumVars(Polysys &sys);
Reduced convertToSquare(Polysys sys, std::vector<double> &sol);

#endif //CPOLYCALC_SQUARE_H
