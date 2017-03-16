#include <iostream>
#include <complex>
#include "print.h"
#include "Reduced.h"
#include "square.h"
#include "wolfram.h"
#include "parser.h"
#include "normalize.h"

using namespace std;

Reduced::Reduced() {
}

string Reduced::toWolfram() {
	ostringstream s;
	s << "{" << this->linears << endl;
	for (auto &p : this->squares) {
		auto &xy = p.first;
		double&k = p.second;
		s << print::valWithSign(k) << "*x" << xy.second << "^2-x"<<xy.first<<","<<endl;
	}
	s << "}";
	return s.str();
}


Point Reduced::updateNewVars(Point sol) {
    map<int, double> vars;
    for (unsigned i = 0; i < sol.size(); i++) {
        vars[i] = sol[i];
    }
    cout << sol << endl;
    while (vars.size()<this->numVars) {
        for (Square sq : squares) {
            if (!vars.count(sq.first.first) && vars.count(sq.first.second)) {
                vars[sq.first.first] = sq.second * vars[sq.first.second] * vars[sq.first.second];
            }
            if (vars.count(sq.first.first) && !vars.count(sq.first.second)) {
                vars[sq.first.second] = sqrt(vars[sq.first.first] / sq.second);
            }
        }
        for (auto&linear : linears) {
            int numMissing = 0;
            int missingVar;
            double s = 0;
            for (auto&kv : linear) {
                if (kv.first == -1) {
                    s += kv.second;
                    continue;
                }
                if (!vars.count(kv.first)) {
                    numMissing++;
                    if (numMissing > 1)break;
                    missingVar = kv.first;
                } else {
                    s += kv.second * vars[kv.first];
                }
            }
            if (numMissing == 1) {
                vars[missingVar] = s / -linear[missingVar];
            }
        }
    }
    Point ret;
    ret.resize(this->numVars);
    for (auto&kv : vars) {
        ret[kv.first] = kv.second;
    }
    return ret;
}

Replacement getOrCreateSquare(map<pair<int, int>, double> &squares, int x, int &varCounter) {
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
    rep.k = 1;
    rep.x = x;
    rep.y = y;
    auto key = pair<int, int>(rep.y, rep.x);
    squares[key] = 1.;
    return rep;
};

Replacement getOrCreateDiff(Polysys &sys, int x, int y, int &varCounter) {
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
    rep.x = x;
    rep.y = y;
    rep.z = z;
    rep.k = 1;
    return rep;
};

Replacement getOrCreatePower(Polysys &sys, map<pair<int, int>, double> &squares, int x,
        int n, int &varCounter) {
    Replacement rep;
    if (n == 2) {
        return getOrCreateSquare(squares, x, varCounter);
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
    Polynome poly;
    auto yKey = multiset<int>();
    yKey.insert(y);
    poly[yKey] = 1;
    poly[genMonomePowerKey(x, n)] = -1.;
    sys.push_back(poly);
    rep.y = y;
    rep.k = 1;
    return rep;
};

bool reducePower(Polynome &poly, Polysys &sys, map<pair<int, int>, double> &squares,
        int &varCounter) {
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
                pwr = getOrCreatePower(sys, squares, v, m, varCounter);
                replacement.clear();
                replacement.insert(pwr.y);
            } else if (powerStats[v] % 2 == 0) { // x^(2*n)
                m = n / 2;
                n = powerStats[v];
                pwr = getOrCreatePower(sys, squares, v, m, varCounter);
                replacement.clear();
                replacement.insert(pwr.y);
                replacement.insert(pwr.y);
            } else if (powerStats[v] >= 2) {
                m = n;
                n = (powerStats[v] - 1);
                pwr = getOrCreatePower(sys, squares, v, m, varCounter);
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

bool reduceProduct(Polynome &poly, Polysys &sys, map<pair<int, int>, double> &squares,
        int &varCounter) {
    for (auto &pair : poly) {
        auto key = pair.first;
        if (key.size() < 2) continue;
        multiset<int>::iterator i = key.begin();
        int x = *i;
        i++;
        int y = *i;
        Replacement pRes = getOrCreateSquare(squares, x, varCounter);
        Replacement qRes = getOrCreateSquare(squares, y, varCounter);
        Replacement tRes = getOrCreateDiff(sys, x, y, varCounter);
        Replacement lRes = getOrCreateSquare(squares, tRes.z, varCounter);
        double k0 = poly[key];
        poly.erase(key);
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

void Reduced::convertToSquare(Polysys sys) {
    map<pair<int, int>, double> &squares = this->squares; // 0,2 -> k     x_0 = k*x_2^2
    int varCounter = calcNumVars(sys);
    // exclude primitive quadratics
    for (size_t p = 0; p < sys.size(); p++) {
        Polynome poly = sys[p];
        Replacement repl;
        bool isq = isPower(poly, 2, repl);
//        cout << "poly=" << poly << ",q=" << isq << endl;
        if (isq) {
            auto key = pair<int, int>(repl.x, repl.y);
            if (squares.count(key) > 0) {
                if (squares[key] != repl.k) {
                    cout << "inconsystent system" << endl;
                    throw "inconsisten system";
                }
            } else {
                squares[key] = repl.k;
            }
            sys.erase(sys.begin() + p);
            p--;
        }
    }
    cout << "sys=" << sys << endl;
    cout << "==============" << endl;
    // reduce high powers
    for (size_t p = 0; p < sys.size(); p++) {
        auto &poly = sys[p];
        while (true) {
            //cout << "reducing poly" << print::polyToString(poly) << endl;
            auto res = reducePower(poly, sys, squares, varCounter);
            if (!res) break;
        }
    }
    for (size_t p = 0; p < sys.size(); p++) {
        auto &poly = sys[p];
        while (true) {
            bool res = reduceProduct(poly, sys, squares, varCounter);
            if (!res) break;
        }
    }
    //cout << "Collecting linears and substs" << endl; 
    for (auto p = 0U; p < sys.size(); p++) {
        Polynome &poly = sys[p];
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
                bool res = findSqSubst(this->squares, v, repl);
                if (!res) {
                    subst.linears[v] = poly[vars];
                } else {
                    subst.squares[repl.x] += poly[vars] * repl.k;
                }
            }
        }
        linears.push_back(linear);
        substs.push_back(subst);
    }
    fillSquaresAndRoots(*this);
    this->numVars = varCounter;
};

void Reduced::normalize() {
    auto&sys = linears;
    for (auto m = 0U; m < sys.size(); m++) {
        for (auto n = m + 1; n < sys.size(); n++) {
            auto coeff = calcK(sys[n], sys[m]);
            for (auto &k : sys[m]) {
                sys[n][k.first] -= coeff * k.second;
                if (fabs(sys[n][k.first]) < 1E-15) {
                    sys[n].erase(k.first);
                }
            }
        }
    }
}

Reduced::Reduced(string str) {
    auto cvt = parseWolfram(str);
    Polysys sys = parser::system(cvt);
    this->convertToSquare(sys);
}

