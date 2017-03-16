//
// Created by steve on 11.11.2015.
//
#include <stdlib.h>
#include <map>
#include <iostream>
#include "parser.h"
#include "util.h"
#include "interval.h"

using namespace std;

Monome parseMonome(string monomeStr) {
    multiset<int> vars;
    auto k = 1.;
    auto parts = split(monomeStr, "*");
    for (auto part:parts) {
        if (part.find("x") == 0) {
            vars.insert(atoi(part.substr(1).c_str()));
        } else {
            k *= atof(part.c_str());
        }
    }
    return Monome(k, vars);
};

string replace(string src, string sep, string glue) {
    return join(split(src, sep), glue);
}

Polynome parsePoly(string polyStr) {
    auto noSpaced = replace(polyStr, " ", "");
    //cout << "polyStr="<<polyStr << endl;
    auto plusSep = replace(noSpaced, "+", " ");
    //cout << "polyStr="<<polyStr << endl;
    auto minusSep = replace(plusSep, "-", " -");
	auto ms1 = replace(minusSep, "-x", "-1*x");

    //cout << "polyStr="<<polyStr << endl;
    //cout << "minusSep="<<minusSep<< endl;
    polyStr = trim(ms1);
    //cout << "parsing poly:"<<polyStr << endl;

    map<multiset<int>, double> ret;
    auto polyParts = split(polyStr, " ");
    for (auto part : polyParts) {
        auto monome = parseMonome(part);
        auto vars = monome.second;
        ret[vars] += monome.first;
    }
    return ret;
};

namespace parser {
    Polysys system(string systemStr) {
        deque<Polynome> ret;
        auto sysParts = split(systemStr, ";");
        for (auto part:sysParts) {
            ret.push_back(parsePoly(part));
        }
        return ret;
    };
}