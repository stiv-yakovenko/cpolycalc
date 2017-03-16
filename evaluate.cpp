//
// Created by steve on 13.11.2015.
//

#include <math.h>
#include <vector>

using namespace std;

#include "evaluate.h"
#include "interval.h"
#include "algebra.h"

namespace eval {

    double sqInt(double k, int x, int y, vector<interval> &iv) {
        return pow((iv[y].first + iv[y].second) * .5, 2) * k
                - .5 * (iv[x].first + iv[x].second);
    };

    int linearInt(Linear &eq, vector<interval> &iv, double &discrepancy) {
        double ret = (eq.count(-1) > 0) ? eq[-1] : 0;
        double largestInt = 0;
        int largestIntIdx = -1;
        for (auto key : eq) {
            int v = key.first;
            if (v == -1) continue;
            if (iv[v].second - iv[v].first > largestInt) {
                largestInt = iv[v].second - iv[v].first;
                largestIntIdx = v;
            }
            ret += eq[v] * .5 * (iv[v].first + iv[v].second);
        }
        discrepancy = largestInt;
        return largestIntIdx;
    };

    double poly(Polynome &poly, Point& sol) {
        double res = 0;
        for (auto key : poly) {
            res += key.second * monome(key.first, sol);
        }
        return res;
    };

    double linear(Linear &lin, Point& sol) {
        double res = 0;
        for (auto &key : lin) {
            if (key.first == -1) {
                res += key.second;
                continue;
            }

            res += key.second * sol[key.first];
        }
        return res;
    };
	double linear(FastLinear &lin, Point& sol) {
		double res = 0;
		for (auto &key : lin) {
			if (key.first == -1) {
				res += key.second;
				continue;
			}

			res += key.second * sol[key.first];
		}
		return res;
	};

    double monome(multiset<int> k, Point & sol) {
        double res = 1;
        for (int i : k) {
            res *= sol[i];
        }
        return res;
    };

}