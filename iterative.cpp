#include <vector>
#include <iostream>

using namespace std;

#include "algebra.h"
#include "interval.h"
#include "iterative.h"
#include "cubic.h"
#include "print.h"

Point2D projectToQuad(Point &sol, int x, int y, double k) {
    auto c = (.5 - k * sol[y]) / k / k;
    vector<double> xClosest = solveCubic(0, c, -.5 * sol[x] / k / k);
    int bestI = -1;
    double bestX, bestY;
    double bestD;
    for (auto i = 0U; i < xClosest.size(); i++) {
        double xx = xClosest[i];
        double d = sqrt(pow(xx - sol[x], 2) + pow(k * xx * xx - sol[y], 2));
        if (bestI == -1 || bestD > d) {
            bestI = i;
            bestX = xx;
            bestY = k * xx * xx;
            bestD = d;
        }
    }
    return Point2D(bestX, bestY);
};

void intersectQuad(Point&sol, Point &p, int x, int y, double k) {
    auto px = p[x];
    auto py = p[y];
    double a = k * px * px;
    double b = 2 * k * sol[x] * px - py;
    double c = k * sol[x] * sol[x] - sol[y];
    interval in = solveSquare(a, b, c);
    double t1 = in.first;
    double t2 = in.second;
    double x1 = sol[x] + t1*px;
    double x2 = sol[x] + t2*px;
    double d1 = fabs(x1 - sol[x]);
    double d2 = fabs(x2 - sol[x]);
    double tBest = d1 < d2 ? t1 : t2;
    for (auto i = 0U; i < sol.size(); i++) {
        sol[i] += tBest * p[i];
    }
}

void projectToLinear(Linear &linear, Point&vect) {
    double nom = 0, n2 = 0;
    for (auto &k : linear) {
        if (k.first == -1) {
            nom += linear[-1];
        } else {
            nom += k.second * vect[k.first];
            n2 += k.second * linear[k.first];
        }
    }
    double t = -nom / n2;
    for (auto&k : linear) {
        if (k.first == -1) continue;
        vect[k.first] += t * linear[k.first];
    }
};

double cosCoeff(Linear &n, Point &v) {
    double n2 = 0., vn = 0.;
    for (auto &k : n) {
        if (k.first == -1) continue;
        vn += k.second * v[k.first];
        n2 += pow(k.second, 2);
    }
    return vn / n2;
}

void iterativeStep(Reduced& reduced, Point& vect) {
    for (int s = reduced.linears.size() - 1; s >= 0; s--) {
        Linear &eq = reduced.linears[s];
        projectToLinear(eq, vect);
    }
    for (auto& square : reduced.squares) {
        int y = square.first.first;
        int x = square.first.second;
        double k = square.second;
        Point2D proj = projectToQuad(vect, x, y, k);
        vect[x] = proj.first;
        vect[y] = proj.second;
    }
}

bool iterCalc(Reduced &reduced, Point &vect, double eps) {
    int iCnt = 0;
    Point prevP = vect;
    while (true) {
        iCnt++;
        iterativeStep(reduced, vect);
        Point dv = vect - prevP;
        auto l = sqrt(length2(dv));
        if (l < eps) {
            return true;
        }
        if (iCnt % 9000 == 0) {
            cout << "cnt=" << iCnt << ", vect = " << print::point(vect) << endl;
        }
    }
};



