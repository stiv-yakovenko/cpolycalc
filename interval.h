#ifndef CPOLYCALC_INTERVAL_H
#define CPOLYCALC_INTERVAL_H

#include <utility>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include "util.h"

using namespace std;

typedef pair<double, double> interval;
interval narrow(interval &old, interval &n);
bool inside(interval &parent, interval &child);
interval join(interval i1, interval i2);
interval sqrtInterval(interval i) ;

interval addInterval(interval i1, interval i2);
interval operator+(interval i1,interval i2);
interval operator-(interval i1,interval i2);

interval subInterval(interval i1, interval i2);

bool inside(double x, interval &i);

interval addConst(interval i, double c);

interval multiply(interval i, double c);
interval square(interval i);

interval intersect(interval i1, interval i2);
interval reduceInterval(interval &i, interval &src);
interval operator+ (interval &other);

bool isEmpty(interval i);
bool equal(interval i1, interval i2) ;
interval quadInt(double a, double b, interval xInt);
vector<interval> quadRoot(double a, double b, interval c, interval xInt);
interval solveSquare(double a, double b, double c);

#endif //CPOLYCALC_INTERVAL_H
