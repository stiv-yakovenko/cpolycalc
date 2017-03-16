#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H
#include "algebra.h"

double findNonzeroRadius (Reduced&r,Point&p);
double findZeroRadius(Reduced&r,Point&p);
Polynome decompose(Polynome&pol,Point&p);
void findRoots(Box b,Reduced&r,map<Point,double>&roots);
Point center(Box&b);
double getBoxR(Box &b);
vector<interval> tuneCombination(Point&center, Reduced&r, Box&box, int tunedVar);
#endif /* NEIGHBORHOOD_H */

