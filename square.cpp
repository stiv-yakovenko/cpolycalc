//
// Created by steve on 11.11.2015.
//
#include <map>
#include <math.h>
#include <iostream>
#include <algorithm>
#include "square.h"
#include "evaluate.h"
#include "parser.h"
#include "print.h"

void checkPoly(Polynome&poly, Point&sol) {
    auto disc = eval::poly(poly, sol);
    if (fabs(disc) > 1E-4) {
        cout << "problem" << endl;
    } else {
        //cout << "checkPoly:" << print::polyToString(poly) << "discr = " << disc << endl;
    }
}

int calcNumVars(Polysys &sys) {
    int ret = 0;
    for (Polynome poly : sys) {
        for (auto &k : poly) {
            for (auto x : k.first) {
                ret = max(ret, x);
            }
        }
    }
    return ret + 1;
};

map<int, int> calcPowers(const multiset<int> &vars) {
    map<int, int> stat;
    for (auto &i : vars) {
        stat[i] += 1;
    }
    return stat;
};

// y = k*x^n?
bool isPower(Polynome &poly, int n, Replacement &rep) {
    if (poly.size() == 2) {
        map<multiset<int>, double>::iterator it = poly.begin();
        multiset<int> k1 = it->first;
        it++;
        multiset<int> k2 = it->first;
        multiset<int> longer = k1.size() > k2.size() ? k1 : k2;
        multiset<int> shorter = k1.size() < k2.size() ? k1 : k2;
        map<int, int> powers = calcPowers(longer);
        if (powers.size() != 1) return false;
        if (shorter.size()!=1) return false;
        if (longer.size() != n) return false;
        rep.x = *shorter.begin();
        rep.y = *longer.begin();
        rep.k = -(poly[longer]) / (poly[shorter]);
        return true;
    }
    return false;
};

Replacement getOrCreateSquare(map<pair<int, int>, double> &squares, int x, int &varCounter, vector<double> &sol) {
    Replacement rep;
    for (Square square : squares) {
        auto vars = square.first;
        if (square.first.second == x) {
            rep.y = vars.first;
            rep.k = square.second;
            return rep;
        }
    }
    auto y = varCounter;
    varCounter++;
    sol.push_back(pow(sol[x], 2));
    rep.k = 1;
    rep.x = x;
    rep.y = y;
    auto key = pair<int, int>(rep.y, rep.x);
    squares[key] = 1.;
    return rep;
};

// z = x-y
Replacement getOrCreateDiff(Polysys &sys, int x, int y, int &varCounter, vector<double> &sol) {
    Replacement rep;
    for (auto p = 0U; p < sys.size(); p++) {
        Polynome poly = sys[p];
        if (poly.size() != 3) continue;
        multiset<int> xIdx;
        multiset<int> yIdx;
        multiset<int> zIdx;
        bool doContinue = false;
        for (auto &i : poly) {
            multiset<int> vars = i.first;
            if (vars.size() != 1) {
                doContinue = true;
                break;
            }
            auto var = *vars.begin();
            if (var == x) xIdx = vars;
            if (var == y) yIdx = vars;
            if (var != x && var != y) zIdx = vars;
        }
        if (doContinue) continue;
        if (xIdx.size() == 0) continue;
        if (yIdx.size() == 0) continue;
        if (zIdx.size() == 0) continue;
        if (xIdx == zIdx || zIdx == yIdx) continue;
        if (poly[xIdx] / poly[yIdx] != -1) continue;
        rep.z = *zIdx.begin();
        rep.k = poly[xIdx] / poly[zIdx];
        return rep;
    }

    int z = varCounter;
    varCounter++;
    sol.push_back(sol[x] - sol[y]);
    Polynome polynome;
    auto k = multiset<int>();
    k.insert(x);
    polynome[k] = 1.;
    k = multiset<int>();
    k.insert(y);
    polynome[k] = -1.;
    k = multiset<int>();
    k.insert(z);
    polynome[k] = -1;
    sys.push_back(polynome);
    checkPoly(polynome, sol);
    rep.x = x;
    rep.y = y;
    rep.z = z;
    rep.k = 1;
    return rep;
};

// x*y -> .5*p+.5*q-.5*l, p=x*x, q=y*y, t=x-y, l=t*t

bool reduceProduct(Polynome &poly, Polysys &sys, map<pair<int, int>, double> &squares,
        int &varCounter, vector<double> &sol) {
    for (auto &pair : poly) {
        auto key = pair.first;
        if (key.size() < 2) continue;
        multiset<int>::iterator i = key.begin();
        int x = *i;
        i++;
        int y = *i;
        Replacement pRes = getOrCreateSquare(squares, x, varCounter, sol);
        Replacement qRes = getOrCreateSquare(squares, y, varCounter, sol);
        Replacement tRes = getOrCreateDiff(sys, x, y, varCounter, sol);
        Replacement lRes = getOrCreateSquare(squares, tRes.z, varCounter, sol);
        double k0 = poly[key];
        poly.erase(key);
        //string str = "Replacing: ";
        stringstream str;
        str << "Replacing: ";
        str << "x" << x << "*x" << y << " -> " << print::valWithSign(.5 * pRes.k) << "*x" << pRes.y;
        str << " " << print::valWithSign(.5 * qRes.k) << "*x" << qRes.y;
        str << " " << print::valWithSign(-.5 * lRes.k * tRes.k * tRes.k) << "*x" << lRes.y << "; ";
        str << "x" << pRes.y << "=" << pRes.k << "x" << x << "^2; ";
        str << "x" << qRes.y << "=" << qRes.k << 'x' << y << "^2; ";
        str << "x" << tRes.z << "=" << print::valWithSign(tRes.k) << "*x" << x;
        str << print::valWithSign(-tRes.k) << "*x" << y << "; ";
        str << "x" << lRes.y << "=" << lRes.k << 'x' << tRes.z << "^2; ";
        //cout << str.str() << endl;
        multiset<int> keyBase = key;
        keyBase.erase(x);
        keyBase.erase(y);
        auto k1 = keyBase;
        auto k2 = keyBase;
        auto k3 = keyBase;
        k1.insert(pRes.y);
        k2.insert(qRes.y);
        k3.insert(lRes.y);
        poly[k1] += .5 * pRes.k * k0;
        poly[k2] += .5 * qRes.k * k0;
        poly[k3] += -.5 * lRes.k * tRes.k * tRes.k * k0;
        if (poly[k1] == 0) poly.erase(k1);
        if (poly[k2] == 0) poly.erase(k2);
        if (poly[k3] == 0) poly.erase(k3);
        return true;
    }
    return false;
};

void replaceVar(Polynome &poly, int v, int y, const multiset<int> &key,
        multiset<int> replacement, double k) {
    auto oldK = poly[key];
    auto k1 = key;
    k1.erase(v);
    for (auto i : replacement) {
        k1.insert(i);
    }
    poly.erase(key);
    poly[k1] = oldK * k;
};

multiset<int> genMonomePowerKey(int x, int n) {
    multiset<int> ret;
    for (int i = 0; i < n; i++) {
        ret.insert(x);
    }
    return ret;
};

// y = k (x_i)^n
Replacement getOrCreatePower(Polysys &sys, map<pair<int, int>, double> &squares, int x,
        int n, int &varCounter, vector<double> &sol) {
    Replacement rep;
    if (n == 2) {
        return getOrCreateSquare(squares, x, varCounter, sol);
    }
    for (auto p = 0U; p < sys.size(); p++) {
        auto poly = sys[p];
        bool pwr = isPower(poly, n, rep);
        if (pwr) {
            return rep;
        }
    }
    int y = varCounter;
    varCounter++;
    sol.push_back(pow(sol[x], n));
    Polynome poly;
    auto yKey = multiset<int>();
    yKey.insert(y);
    poly[yKey] = 1;
    poly[genMonomePowerKey(x, n)] = -1.;
    checkPoly(poly, sol);
    sys.push_back(poly);
    rep.y = y;
    rep.k = 1;
    return rep;
};

bool reducePower(Polynome &poly, Polysys &sys, map<pair<int, int>, double> &squares,
        int &varCounter,
        vector<double> &sol) {
    // vars = [5,5,5,8,8,9] = x_5^3 + x_8^2 + x_9
    for (auto &key : poly) {
        if (key.first.size() < 2) continue;
        map<int, int> powerStats = calcPowers(key.first);
        for (auto &v : key.first) {
            int m, n = 0;
            Replacement pwr;
            multiset<int> replacement;
            if (powerStats[v] == 2) {
                m = 2;
                pwr = getOrCreatePower(sys, squares, v, m, varCounter, sol);
                replacement.clear();
                replacement.insert(pwr.y);
            } else if (powerStats[v] % 2 == 0) { // x^(2*n)
                m = n / 2;
                n = powerStats[v];
                pwr = getOrCreatePower(sys, squares, v, m, varCounter, sol);
                replacement.clear();
                replacement.insert(pwr.y);
                replacement.insert(pwr.y);
            } else if (powerStats[v] >= 2) {
                m = n;
                n = (powerStats[v] - 1);
                pwr = getOrCreatePower(sys, squares, v, m, varCounter, sol);
                replacement.clear();
                replacement.insert(pwr.y);
                replacement.insert(v);
            } else {
                continue;
            }
            //cout << "replace x" << pwr.y << "=" << pwr.k << "*x" << v << "^" << m << endl;
            int y = pwr.y;
            replaceVar(poly, v, y, key.first, replacement, pwr.k);
            return true;
        }
    }
    return false;
};

bool findSqSubst(std::map<std::pair<int, int>, double> &squares, int y, Replacement &repl) {
    for (auto key : squares) {
        if (key.first.first == y) {
            repl.x = key.first.second;
            repl.k = key.second;
            return true;
        }
    }
    return false;
};

void fillSquaresAndRoots(Reduced&reduced) {
    for (auto &square : reduced.squares) {
        XY xy = square.first;
        int x = xy.second;
        int y = xy.first;
        double k = square.second;
        reduced.xToKY[x] = std::pair<int, double>(y, k);
    }
}

Reduced convertToSquare(Polysys sys, vector<double> &sol) {
    Reduced reduced; // = {squares: {}, linears: sys, substs: []};
    map<pair<int, int>, double> &squares = reduced.squares; // 0,2 -> k     x_0 = k*x_2^2
    int varCounter = calcNumVars(sys);
    // exclude primitive quadratics
    for (size_t p = 0; p < sys.size(); p++) {
        Polynome poly = sys[p];
        Replacement repl;
        bool isq = isPower(poly, 2, repl);
        //cout<<"poly="<<poly<<",isq="<<isq<<endl;
        if (isq) {
            auto key = pair<int, int>(repl.x, repl.y);
            if (squares.count(key) > 0) {
                if (squares[key] != repl.k) {
                    throw "inconsisten system";
                }
            } else {
                squares[key] = repl.k;
            }
            sys.erase(sys.begin() + p);
            p--;
        }
    }
    cout <<"sys"<<sys<<endl;
    // reduce high powers
    for (size_t p = 0; p < sys.size(); p++) {
        auto &poly = sys[p];
        while (true) {
            //cout << "reducing poly" << print::polyToString(poly) << endl;
            auto res = reducePower(poly, sys, squares, varCounter, sol);
            if (!res) break;
        }
    }
    for (size_t p = 0; p < sys.size(); p++) {
        auto &poly = sys[p];
        while (true) {
            bool res = reduceProduct(poly, sys, squares, varCounter, sol);
            if (!res) break;
        }
    }
    //cout << "Collecting linears and substs" << endl; 
    for (auto p = 0U; p < sys.size(); p++) {
        Polynome &poly = sys[p];
        checkPoly(poly, sol);
        Subst subst;
        Linear linear;
        for (auto &key : poly) {
            const multiset<int> &vars = key.first;
            if (vars.size() == 0) {
                subst.linears[-1] = poly[vars];
                linear[-1] = poly[vars];
            } else {
                int v = *key.first.begin();
                linear[v] = poly[vars];
                Replacement repl;
                bool res = findSqSubst(reduced.squares, v, repl);
                if (!res) {
                    subst.linears[v] = poly[vars];
                } else {
                    subst.squares[repl.x] += poly[vars] * repl.k;
                }
            }
        }
        reduced.linears.push_back(linear);
        reduced.substs.push_back(subst);
    }

    fillSquaresAndRoots(reduced);

    reduced.numVars = varCounter;

    return reduced;
};





