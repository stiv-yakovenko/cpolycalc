//
// Created by steve on 13.11.2015.
//

#ifndef CPOLYCALC_EVALUATE_H
#define CPOLYCALC_EVALUATE_H

#include <vector>

#include "interval.h"
#include "algebra.h"

namespace eval {
    int linearInt(Linear &eq, std::vector<interval> &iv, double &discrepancy);
    double sqInt(double k, int x, int y, std::vector<interval> &iv);
    double poly(Polynome &poly, Point& sol);
    double monome(multiset<int> k, Point &sol);
    double linear(Linear &lin, Point& sol);
	double linear(FastLinear &lin, Point& sol);
}

#endif //CPOLYCALC_EVALUATE_H
