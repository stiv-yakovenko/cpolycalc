#ifndef REDUCED_H
#define REDUCED_H

#include "algebra.h"
using namespace std;

struct Reduced {
    Reduced();
    Reduced(string sys);
    void convertToSquare(Polysys sys);
    void normalize();
    Point updateNewVars(Point sol);
	string toWolfram();
    
    int numVars;		
    map<int, Polynome> gradSteps;
    Polynome grandPoly;
    map<multiset<int>,Polynome> grandDerivatives;
    deque<Linear> linears;
    map<XY, double> squares;
    map<int, pair<int,double> > xToKY;
    map<int, Linear> quadToLin;
    map<int, Linear> linToQuad;
    set<int> linVars, quadVars;
    
    deque<Subst> substs;
};

#endif /* REDUCED_H */

