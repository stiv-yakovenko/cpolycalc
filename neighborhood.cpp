#include <cstdlib>
#include <cmath>
#include <math.h>
#include <iostream>
#include "print.h"
#include "float.h"
#include "quartic.h"
#include "gradient.h"
#include "neighborhood.h"
#include "interval.h"
#include "evaluate.h"
#include "eigen.h"
#include "calc.h"
#include "poly34.h"
#include "geom.h"

Polynome decompose(Polynome&pol, Point&p) {
    Polynome res;
    for (auto&kv : pol) {
        Polynome product;
        multiset<int> one;
        product[one] = kv.second;
        for (auto&x : kv.first) {
            Polynome mult;
            multiset<int> x_i;
            x_i.insert(x);
            mult[x_i] = 1;
            multiset<int> b;
            mult[b] = p[x];
            product = product*mult;
        }
        res += product;
    }
    return res;
}

double getBoxR(Box &b) {
    double sum = 0.;
    for (auto&pair : b) {
        sum += (pair.second - pair.first)*(pair.second - pair.first);
    }
    return 0.5 * sqrt(sum);
}

Point center(Box&b) {
    Point r;
    for (auto&pair : b) {
        r.push_back((pair.second + pair.first)*.5);
    }
    return r;
}

bool boxInsideRoot(Box&b, Point&p, double&R) {
    double r = getBoxR(b);
    Point c = center(b);
    Point d = p - c;
    return length(d) + r<R;
}

vector<interval> tuneCombination(Point&center, Reduced&r, Box&box, int tunedVar) {
    /*    Point sol__ = Point({1.7634461203153509072715803028927950560639240811099,
            1.3541284792039218991442153033803225083293561086555,
            2.7247999726431325494720518840050925341568038468,
            1.0788127139187973220560762548608949507249579816996,
            -0.28626298910698669344619665590600085824200004863,
            -0.6703209217045077554834247470797184959849024162400});
        Point sol = r.updateNewVars(sol__);
        cout<<"inside="<<inside(box,sol)<<endl;*/
    map<int, Polynome> grads;
    int n = 0;
    for (auto&linear : r.linears) {
        Polynome grad;
        for (auto&kv : linear) {
            if (kv.first == -1) continue;
            Polynome dp;
            multiset<int> x;
            x.insert(n);
            dp[x] = kv.second;
            grads[kv.first] += dp;
        }
        n++;
    }
    multiset<int> key;
    for (auto&sq : r.squares) {
        auto &xy = sq.first;
        auto&k = sq.second;
        Polynome dg;
        key = multiset<int>();
        key.insert(n);
        dg[key] = -2. * k * center[xy.second];
        grads[xy.second] += dg;
        dg = Polynome();
        key = multiset<int>();
        key.insert(n);
        dg[key] = 1.;
        grads[xy.first] += dg;
        n++;
    }
    Polynome one;
    key = multiset<int>();
    one[key] = -1;
    grads[tunedVar] += one;
    Polynome G;
    for (auto&kv : grads) {
        Polynome sq = kv.second * kv.second;
        G += sq;
    }
	//cout << "G=" << G << endl;
    Point p1 = minimizeQForm(G, center);
	for (auto&pp : p1) {
		if (!isfinite(pp)) pp = 1.;
	}
	//cout << "coeffs" << p1 << endl;
    //p1[1]= 
    double A = 0, B = 0;
    Polynome C;
    n = 0;
    for (auto&linear : r.linears) {
        double &coef = p1[n];
        for (auto&kv : linear) {
            if (kv.first == -1) {
                double pr = kv.second*coef;
                C += pr;
            } else if (kv.first == tunedVar) {
                B += kv.second*coef;
            } else {
                Monome m;
                m.second.insert(kv.first);
                m.first = kv.second*coef;
                C += m;
            }
        }
        n++;
    }
    // y-k*x*x
    for (auto&sq : r.squares) {
        double &coef = p1[n];
        auto &yx = sq.first;
        auto&k = sq.second;
        if (yx.second == tunedVar) {
            Square z = sq;
            A += -coef*k;
        } else {
            Monome m;
            m.first = -k*coef;
            m.second.insert(yx.second);
            m.second.insert(yx.second);
            C += m;
        }
        if (yx.first == tunedVar) {
            B += coef;
        } else {
            Monome m;
            m.first = coef;
            m.second.insert(yx.first);
            C += m;
        }
        n++;
    }
    //cout << "A=" << A << endl;
    //cout << "B=" << B << endl;
    //cout << "C=" << C << endl;
    //cout << "C[sol]=" << eval::poly(C, sol) << endl;
    //    cout << "r=" << findNonzeroRadius(r, sol) << endl;
    multiset<int> k_0;
    interval cInt = interval(C[k_0], C[k_0]);
    for (int i = 0; i < r.numVars; i++) {
        if (i == tunedVar)continue;
        multiset<int> k_i;
        k_i.insert(i);
        double b = C[k_i];
        k_i.insert(i);
        double a = C[k_i];
        interval r = quadInt(a, b, box[i]);
        cInt = cInt + r;
    }
    //cout << "cInt=" << cInt << endl;
    //cout<<"A="<<A<<", B="<<B<<endl;
    vector<interval> xRoots = quadRoot(A, B, cInt, box[tunedVar]);
    //cout << "x" << tunedVar << "=" << box[tunedVar] << endl;
    //cout << "xRoots=" << xRoots << endl;
    //    cout<<"inside="<<inside(box,sol)<<endl;
    return xRoots;
}

int numIt = 0;

void findRoots(Box b, Reduced&r, map<Point, double>&roots) {
    cout << "IT:" << numIt << ",roots=" << roots.size() << ", b:" << b << endl;
    for (auto&kv : roots) {
        Point root = kv.first;
        if (boxInsideRoot(b, root, kv.second)) {
            return;
        }
    }
    doTune(r, b);
    numIt++;
    auto R = getBoxR(b);
    Point p;
    int longest = -1, i = 0;
    if (R < 1E-6) {
        cout << "SMALL" << endl;
    }
    for (auto&pair : b) {
        auto d = pair.second - pair.first;
        p.push_back((pair.second + pair.first)*.5);
        if (longest < 0 || d > (b[longest].second - b[longest].first)) {
            longest = i;
        }
        i++;
    }
    auto val = eval::poly(r.grandPoly, p);
    if (val < 1E-4) {
        Point p1 = p;
        descentQuad(r, p1);
        auto val = eval::poly(r.grandPoly, p1);
        if (val < 1E-10) {
            cout << "ROOT" << p << endl;
            auto Rz = findZeroRadius(r, p1);
            roots[p1] = Rz;
        }
    }
    auto nonZeroR = findNonzeroRadius(r, p);
    //    cout << "boxR="<<R<<", nonZeroR="<<nonZeroR<<endl;
    if (R < nonZeroR) {
        return;
    }
    Box left = b;
    Box right = b;
    double mid = .5 * (b[longest].first + b[longest].second);
    left[longest] = interval(b[longest].first, mid);
    right[longest] = interval(mid, b[longest].second);
    findRoots(left, r, roots);
    findRoots(right, r, roots);
}

double findZeroRadius(Reduced&r, Point&p) {
    Polynome t = tailor(r, p);
    int n = r.numVars;
    vector<double> mat(n*n, 0.);
    double A = 0, B = 0, C = 0;
    // A rr + Brrr + Crrrr
    for (auto&kv : t) {
        auto pwr = kv.first.size();
        if (pwr == 2) {
            auto it = kv.first.begin();
            auto i = *it;
            it++;
            auto j = *it;
            double k = i == j ? 1. : .5;
            mat[i * n + j] = kv.second*k;
            mat[j * n + i] = kv.second*k;
        }
        if (pwr == 3) {
            B += -fabs(kv.second);
        }
        if (pwr == 4) {
            C += fabs(kv.second);
        }
    }
    auto vals = eigenValues(mat, n);
    A = *min_element(&vals[0], &vals[0] + n);
    interval res = solveSquare(C, B, A);
    if (res.second < res.first) return 0;
    if (res.first > 0) {
        return res.first;
    }
    return res.second;
}

double findNonzeroRadius(Reduced&r, Point&p) {
    double R = DBL_MAX;
    for (auto& kv : r.squares) {
        auto&xy = kv.first;
        double&k = kv.second;
        R = min(R, distToQuad(p, xy.first, xy.second, k));
    }
    for (auto&l : r.linears) {
        R = min(R, dist(l, p));
    }
    return R;
}


