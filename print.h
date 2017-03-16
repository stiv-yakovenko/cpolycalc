//
// Created by steve on 12.11.2015.
//

#ifndef CPOLYCALC_PRINT_H
#define CPOLYCALC_PRINT_H

#include <string>

#include "Reduced.h"
#include "interval.h"
#include "algebra.h"

using namespace std;

namespace print {
    string monome(const multiset<int> &vars);
    string valWithSign(double v);
    string sys(deque<Linear> &sys, vector<double> &sol);
    string sys(deque<Polynome> &sys, vector<double> &sol);
    string polyToString(Polynome &poly);
    string linToString(Linear&lin, Point&sol);
    string reduced(Reduced&red, Point& sol);
    string point(Point&sol);
    string vect(Point&sol);
    string normal(Linear&linear, Point& sol);
};

ostream& operator<<(ostream& stream, interval& iv);
ostream& operator<<(ostream& stream, Polysys& iv);
ostream& operator<<(ostream& stream, Box& iv);
ostream& operator<<(ostream& stream, Linear& lin);
ostream& operator<<(ostream& stream, Square& sq);
ostream& operator<<(ostream& stream, Point& iv);
ostream& operator<<(ostream& stream, deque<Linear>& iv);
ostream& operator<<(ostream& stream, map<XY, double>& iv);
ostream& operator<<(ostream& stream, const Point& iv);
ostream& operator<<(ostream& stream, Polynome& iv);
ostream& operator<<(ostream& stream, Subst& iv);
ostream& operator<<(ostream& stream, Monome& mon);
ostream& operator<<(ostream& stream, multiset<int>& mon);
ostream& operator<<(ostream& stream, const multiset<int>& mon);
ostream& operator<<(ostream& stream, vector<Point>& points);

#endif //CPOLYCALC_PRINT_H
