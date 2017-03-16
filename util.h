//
// Created by steve on 11.11.2015.
//

#ifndef CPOLYCALC_UTIL_H
#define CPOLYCALC_UTIL_H

#include <vector>
#include <string>
#include <sstream>

using namespace std;

template<typename T>
vector<T> concat(vector<T> &a, vector<T> &b) {
    vector<T> ret = vector<T>();
    copy(a.begin(), a.end(), back_inserter(ret));
    copy(b.begin(), b.end(), back_inserter(ret));
    return ret;
}

vector<string> split(string str, string sep);
string join(vector<string> arr, string glue);
string trim(string &str);

double minElem(vector<double>&vect);
double maxElem(vector<double>&vect);

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif //CPOLYCALC_UTIL_H
