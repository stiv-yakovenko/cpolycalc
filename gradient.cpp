#include <set>
#include <functional>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <float.h>

#include "iterative.h"
#include "evaluate.h"
#include "algebra.h"
#include "print.h"
#include "quartic.h"
#include "evaluate.h"
#include "geom.h"
#include "parser.h"
#include "gradient.h"
#include "cubic.h"
using namespace std;

int getOnePower(const multiset<int> &k, int y) {
    int power = 0;
    for (int p : k) {
        if (p == y) power++;
    }
    return power;
}

void calcGradSteps(Reduced&ren) {
    for (int y = 0; y < ren.numVars; y++) {
        int maxPower = 0;
        map<int, Polynome> coeffs;
        for (auto key : ren.grandPoly) {
            int power = getOnePower(key.first, y);
            maxPower = max(maxPower, power);
            // divide by y^power
            multiset<int> k = key.first;
            for (int p = 0; p < power; p++) {
                k.erase(y);
            }
            if (coeffs.count(power) == 0) {
                coeffs[power] = Polynome();
            }
            Polynome &coeff = coeffs[power];
            coeff[k] += key.second;
        }
        if (maxPower != 2) continue;
        Polynome &a = coeffs[2];
        multiset<int> emptyKey;
        if (a.size() == 1 && a.count(emptyKey) == 1) {
            double A = a[emptyKey];
            ren.gradSteps[y] = coeffs[1]; // B
            for (auto &key : ren.gradSteps[y]) {
                key.second *= -.5 / A;
            }
        }
    }
}

void minimize(Reduced &reduced, Point&sol) {
    int step = 0;
    while (true) {
        step++;
        for (auto &gradKey : reduced.gradSteps) {
            auto y = gradKey.first;
            Polynome &p = gradKey.second;
            double y1 = eval::poly(p, sol);
            sol[y] = y1;
        }
        for (auto &square : reduced.squares) {
            double k = square.second;
            int y = square.first.first;
            int x = square.first.second;
            double y0 = sol[y] / k;
            if (y0 < 0) {
                sol[x] /= 2;
                sol[y] = k * sol[x] * sol[x];
                continue;
            }
            double x1 = -sqrt(y0);
            double x2 = +sqrt(y0);
            sol[x] = x1;
            double min1 = eval::poly(reduced.grandPoly, sol);
            sol[x] = x2;
            double min2 = eval::poly(reduced.grandPoly, sol);
            if (min1 < min2) {
                sol[x] = x1;
            } else {
                sol[x] = x2;
            }
        }
        if (step % 100 == 0) {
            double g = eval::poly(reduced.grandPoly, sol);
            cout << endl << endl << endl;
            cout << "g=" << g << endl;
            cout << print::point(sol) << endl;
        }
    }
}

void calcG(Reduced &r) {
    Polynome &p = r.grandPoly;
    for (Linear &linear : r.linears) {
        for (auto &k1 : linear) {
            for (auto &k2 : linear) {
                multiset<int> k;
                if (k1.first != -1) k.insert(k1.first);
                if (k2.first != -1) k.insert(k2.first);
                double v = k1.second * k2.second;
                p[k] += v;
            }
        }
    }

    for (Square sq : r.squares) {
        int y = sq.first.first;
        int x = sq.first.second;
        double k = sq.second;
        // (y-k*x*x)^2 = y^2 - 2k*y*x^2 + k*k*x*x*x*x
        multiset<int> k1;
        k1.insert(y);
        k1.insert(y);
        p[k1] += 1;
        multiset<int> k2;
        k2.insert(y);
        k2.insert(x);
        k2.insert(x);
        p[k2] += -2 * k;
        multiset<int> k3;
        k3.insert(x);
        k3.insert(x);
        k3.insert(x);
        k3.insert(x);
        p[k3] += k*k;
    }
}

Monome monomeDeriv(const multiset<int> &num, multiset<int>&denum) {
    map<int, int> powers;
    for (auto&p : num) {
        if (powers.count(p)) {
            powers[p]++;
        } else {
            powers[p] = 1;
        }
    }
    double k = 1.;
    auto ret = Monome(0., multiset<int>());
    for (auto&p : denum) {
        if (!powers.count(p)) {
            return ret;
            break;
        }
        k *= powers[p];
        powers[p]--;
    }
    for (auto&kv : powers) {
        for (int i = 0; i < kv.second; i++) {
            ret.second.insert(kv.first);
        }
    }
    ret.first = k;
    return ret;
}

set<int> getAllVars(Reduced&r) {
    set<int> vars;
    for (auto &kv : r.grandPoly) {
        const multiset<int> &mvars = kv.first;
        for (auto &v : mvars) {
            vars.insert(v);
        }
    }
    return vars;
}

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void iterateQuad(set<int> &vars, function<void(multiset<int>) > cb) {
    for (auto &i0 : vars) for (auto &i1 : vars) {
            multiset<int> k1;
            k1.insert(i0);
            k1.insert(i1);
            cb(k1);
        }
    for (auto &i0 : vars) {
        multiset<int> k1;
        k1.insert(i0);
        cb(k1);
    }
    cb(multiset<int>());
}

void iterateFast(multiset<int>usedVars, Polynome p, Reduced&r) {
    if (usedVars.size()>=4)return;
    for (int i = -1; i < r.numVars; i++) {
        multiset<int> k1;
        if (i>-1) k1.insert(i);        
        multiset<int> vars(usedVars);
        if (i>-1) vars.insert(i);
        if (r.grandDerivatives.count(vars)>0) continue;
        Polynome &der = r.grandDerivatives[vars];
        for (auto &kv :p) {
            Monome monDer = monomeDeriv(kv.first, k1);
            multiset<int> p = kv.first;
            if (monDer.first == 0) continue;
            monDer.first *= kv.second;
            add(der, monDer);
        }        
        removeZeroes(der);
        iterateFast(vars,der,r);
    }
}
void calcDerivativesFast(Reduced&r){
    iterateFast(multiset<int>(),r.grandPoly,r);
}
// ugly buy fast-to-implement
void iterate(set<int> &vars, function<void(multiset<int>) > cb) {
    for (auto &i0 : vars) for (auto &i1 : vars) for (auto &i2 : vars) for (auto &i3 : vars) {
                    multiset<int> k1;
                    k1.insert(i0);
                    k1.insert(i1);
                    k1.insert(i2);
                    k1.insert(i3);
                    cb(k1);
                }
    for (auto &i0 : vars) for (auto &i1 : vars) for (auto &i2 : vars) {
                multiset<int> k1;
                k1.insert(i0);
                k1.insert(i1);
                k1.insert(i2);
                cb(k1);
            }
    for (auto &i0 : vars) for (auto &i1 : vars) {
            multiset<int> k1;
            k1.insert(i0);
            k1.insert(i1);
            cb(k1);
        }
    for (auto &i0 : vars) {
        multiset<int> k1;
        k1.insert(i0);
        cb(k1);
    }
    cb(multiset<int>());
}

Polynome tailor(Reduced&r, Point&p) {
    auto vars = getAllVars(r);
    Polynome ret;
    iterate(vars, [&r, &p, &ret](multiset<int>k1)->void {
        Polynome &der = r.grandDerivatives[k1];
        double a = eval::poly(r.grandDerivatives[k1], p) / factorial(k1.size());
        if (a == 0.)return;
            if (ret.count(k1)) {
                ret[k1] += a;
            } else {
                ret[k1] = a;
            }
    });
    return ret;
}

void calcGDerivatives(Reduced &r) {
    auto vars = getAllVars(r);
    cout << "calculating derivatives" << endl;
    iterate(vars, [&r](multiset<int>k1)->void {
        if (r.grandDerivatives.count(k1)) return;
        Polynome der = r.grandDerivatives[k1];
        for (auto &kv : r.grandPoly) {
            Monome monDer = monomeDeriv(kv.first, k1);
                    multiset<int> p = kv.first;
            if (monDer.first == 0) continue;
                    monDer.first *= kv.second;
                    add(der, monDer);
            }
        removeZeroes(der);
                r.grandDerivatives[k1] = der;
    });
}

Polynome dirPoly(Reduced&r, Point&p, Point&v) {
    Polynome res;
    for (auto &kv : r.grandPoly) {
        Polynome product;
        product[multiset<int>()] = kv.second;
        for (auto&x : kv.first) {
            Polynome pvt; // p_x + v_x*t, t==x_0
            Monome p_x;
            p_x.first = p[x];
            Monome v_x_t;
            v_x_t.second.insert(0);
            v_x_t.first = v[x];
            add(pvt, p_x);
            add(pvt, v_x_t);
            product = mult(product, pvt);
        }
        res = res + product;
    }
    removeZeroes(res);
    return res;
}

void optimize(Point &v, Reduced&r, Point&p) {
    Polynome pol = dirPoly(r, p, v);
    double c[5];
    for (auto&kv : pol) {
        c[kv.first.size()] = kv.second;
    }
/*    q[0] = 3. * c[3] / c[4] / 4.;
    q[1] = 2. * c[2] / c[4] / 4.;
    q[2] = c[1] / c[4] / 4.;*/
    auto roots = solveCubic(c[1] / c[4] / 4.,2. * c[2] / c[4] / 4.,3. * c[3] / c[4] / 4.);
    int n = roots.size();
    double bestVal = DBL_MAX;
    int bestI = -1;
    for (int i = 0; i < n; i++) {
        Point t;
        t.push_back(roots[i]);
        auto val = eval::poly(pol, t);
        if (val < bestVal || bestI == -1) {
            bestVal = val;
            bestI = i;
        }
    }
    if (bestI == -1) return;
    auto t = roots[bestI];
    //    cout << "t="<<t<<endl;
    Point v1 = v *t;
    p += v1;
}

int cnt = 0;

Point grad(Reduced&r, Point&p) {
    Point v;
    for (unsigned i = 0; i < p.size(); i++) {
        multiset<int> k;
        k.insert(i);
        Polynome &gr = r.grandDerivatives[k];
        v.push_back(eval::poly(gr, p));
    }
    return v;
}

void descent(Reduced&r, Point &p0) {
    cout << "Gradient descent" << endl;
    Point p = p0;
    int n = getAllVars(r).size();
    auto val = eval::poly(r.grandPoly, p);
    while (true) {
        cnt++;
        if (cnt % 1000 == 0) {
            cout << "CNT:" << cnt << "val=" << val << ", v=" << p << endl;
        }
        Point v = grad(r, p);
        v=v/length(v);
        optimize(v, r, p);
        auto newVal = eval::poly(r.grandPoly, p);
        if (newVal >= val) {
            cout << "Finished optimization" << endl;
            return;
        }
        val = newVal;
    }
}
void descentQuad(Reduced&r, Point &p) {
    cout << "Gradient descent" << endl;
    int n = getAllVars(r).size();
    auto val = eval::poly(r.grandPoly, p);
    while (true) {
        cnt++;
        cout << "CNT:" << cnt << "G=" << val << ", v=" << p << endl;
        Point v = grad(r, p);
        double l = length(v);
        if(l!=0.){
            v=v/l;
            optimize(v, r, p);
        }
        Polynome t=quadTailor(r, p);
        Point v1 = minimizeQForm(t, p);
        optimize(v1, r, p);
        auto newVal = eval::poly(r.grandPoly, p);
        if (newVal >= val) {
            cout << "Finished optimization" << endl;
            return;
        }
        val = newVal;
    }
}

Point sqProject(Reduced& reduced, Point vect) {
    for (auto& square : reduced.squares) {
        int y = square.first.first;
        int x = square.first.second;
        double k = square.second;
        Point2D proj = projectToQuad(vect, x, y, k);
        vect[x] = proj.first;
        vect[y] = proj.second;
    }
    return vect;
}

Point lineProject(Reduced& reduced, Point vect) {
    for (int s = reduced.linears.size() - 1; s >= 0; s--) {
        Linear &eq = reduced.linears[s];
        projectToLinear(eq, vect);
    }
    return vect;
}

Polynome quadTailor(Reduced&r, Point &p) {
    auto vars = getAllVars(r);
    Polynome ret;
    iterateQuad(vars, [&r, &p, &ret](multiset<int>k1)->void {
        Polynome &der = r.grandDerivatives[k1];
        double a = eval::poly(r.grandDerivatives[k1], p) / factorial(k1.size());
        if (a == 0.)return;
            if (ret.count(k1)) {
                ret[k1] += a;
            } else {
                ret[k1] = a;
            }
        if (ret[k1] == 0.) ret.erase(k1);
        });
    return ret;
}

Point minimizeQForm(Polynome qTailor, Point p) {
    Point p0 = p;
    vector<pair<int, Polynome>> eqs;
    for (unsigned i = 0; i < p.size(); i++) {
        ////        cout<<"tailor="<<qTailor<<endl;
        //     cout << "excluding x_"<<i <<endl;
        // ax_i^2 + bx_i = a(x_i+b/(2a))^2 - b*b/4/a/a
        // R = - b*b/4/a
        double a = 0;
        Polynome b;
        set<multiset<int>> toRemove;
        for (auto&kv : qTailor) {
            //       cout << "term " << kv.first<<endl;
            auto cnt = kv.first.count(i);
            if (cnt == 2) {
                a = kv.second;
                //         cout<<"kv.second="<<kv.second<<endl;
                //       cout<<"a="<<a<<endl;
                toRemove.insert(kv.first);
            } else if (cnt == 1) {
                auto b_i = kv.first;
                b_i.erase(i);
                b[b_i] += kv.second;
                toRemove.insert(kv.first);
            }
        }
        //cout << "b="<<b<<",a="<<a<<endl;
        for (auto&k : toRemove) {
            qTailor.erase(k);
        }
        Polynome R = b*b;
        R *= .25 / a;
        //cout<<"R="<<R<<endl;
        qTailor -= R;
        b *= -.5 / a;
        eqs.push_back(pair<int, Polynome>(i, b));
        //cout<<"x"<<i<<"="<<b<<endl;        
    }
    for (auto i = eqs.rbegin(); i != eqs.rend(); ++i) {
        p[i->first] = eval::poly(i->second, p);
    }
    return p;
}
/*
int main(){
    Polynome pol=parsePoly("211*x0*x0+x1*x1+3*x0*x1+2*x0-33*x1");
    Point p=Point({0,0});
    minimizeQForm(pol,p);
    cout<<p<<endl;    
}
 */