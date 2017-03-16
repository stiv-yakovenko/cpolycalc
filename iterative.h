#ifndef ITERATIVE_H
#define ITERATIVE_H


#include "algebra.h"
#include "2d.h"
#include "Reduced.h"

void iterativeStep(Reduced& reduced, Point& vect);
bool iterCalc(Reduced &reduced, Point &vect,double eps);
void projectToLinear(Linear &linear, Point&vect);
Point2D projectToQuad(Point &sol, int x, int y, double k);
void intersectQuad(Point&sol, Point &p, int x, int y, double k);

#endif