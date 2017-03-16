//
// Created by steve on 13.11.2015.
//

#ifndef CPOLYCALC_CALC_H
#define CPOLYCALC_CALC_H

#include "algebra.h"
#include "interval.h"
#include "Reduced.h"
int doTune(Reduced &reduced, Box &box);
double calcV(vector<interval> &iv);
bool boxCalc(Reduced &reduced, std::vector<interval> iv, int depth, std::vector<double> &sol, std::vector<Point> &roots, double progress);

#endif //CPOLYCALC_CALC_H
