#include "math.h"
#include "algebra.h"

double operator*(Point&p1, Point&p2) {
    double d = 0;
    for (auto i = 0U; i < p1.size(); i++) {
        d += p1[i] * p2[i];
    }
    return d;
}

Point operator*(Point&p1, double c) {
    Point p = p1;
    for (auto i = 0U; i < p1.size(); i++) {
        p[i] *= c;
    }
    return p;
}

void operator+=(Point&me, Point&d) {
    for (auto i = 0U; i < me.size(); i++) {
        me[i] += d[i];
    }
}

void operator-=(Point&me, Point&d) {
    for (auto i = 0U; i < me.size(); i++) {
        me[i] -= d[i];
    }
}

Point operator/(Point&p1, double c) {
    Point p = p1;
    for (auto i = 0U; i < p1.size(); i++) {
        p[i] /= c;
    }
    return p;
}

Point operator+(Point&p1, Point&p2) {
    Point p = p1;
    for (auto i = 0U; i < p1.size(); i++) {
        p[i] += p2[i];
    }
    return p;
}

void subtract(Point &x1, Point &x2) {
    for (unsigned int i = 0; i < x1.size(); i++) {
        x1[i] -= x2[i];
    }
}

double length2(Point&p) {
    double d = 0;
    for (unsigned int i = 0; i < p.size(); i++) {
        d += p[i] * p[i];
    }
    return d;
}

Point sub(Point&p1, Point&p2) {
    Point r;
    for (unsigned int i = 0; i < p1.size(); i++) {
        r.push_back(p1[i] - p2[i]);
    }
    return r;
}

Point operator-(Point&p1, Point&p2) {
    return sub(p1, p2);
}

double length(Point&p) {
    return sqrt(length2(p));
}

bool inside(Box&box, Point&p) {
    for (auto i = 0U; i < p.size(); i++) {
        if (p[i] < box[i].first || box[i].second < p[i]) return false;
    }
    return true;
}

void add(Polynome&poly, Monome&mon) {
    if (poly.count(mon.second)) {
        poly[mon.second] += mon.first;
    } else {
        poly[mon.second] = mon.first;
    }
}

Monome mult(Monome&a, Monome&b) {
    Monome ret;
    ret.first = a.first * b.first;
    ret.second.insert(a.second.begin(), a.second.end());
    ret.second.insert(b.second.begin(), b.second.end());
    return ret;
}

void operator/=(Polynome&p, double c) {
    for (auto &kv : p) {
        kv.second /= c;
    }
}

void operator-=(Polynome&p, Polynome&s) {
    for (auto&kv : s) {
        p[kv.first] -= kv.second;
        if (p[kv.first] == .0) p.erase(kv.first);
    }
}

void operator+=(Polynome&p, Polynome&s) {
    for (auto&kv : s) {
        p[kv.first] += kv.second;
        if (p[kv.first] == .0) p.erase(kv.first);
    }
}

void operator*=(Polynome&p, double c) {
    for (auto &kv : p) {
        kv.second *= c;
    }
}

Polynome mult(Polynome&a, Polynome&b) {
    Polynome ret;
    for (auto &kv : a) {
        for (auto &kv1 : b) {
            Monome m1 = Monome(kv.second, kv.first);
            Monome m2 = Monome(kv1.second, kv1.first);
            Monome prod = mult(m1, m2);
            ret[prod.second] += prod.first;
        }
    }
    return ret;
}

Polynome operator*(Polynome&a, Polynome&b) {
    return mult(a, b);
}

void removeZeroes(Polynome&p) {
    set<multiset<int>> zeroes;
    for (auto &kv : p) {
        if (p[kv.first] == 0.) {
            zeroes.insert(kv.first);
        }
    }
    for (auto&k : zeroes) {
        p.erase(k);
    }
}

Polynome operator+(Polynome&a, Polynome&b) {
    Polynome ret;
    for (auto &kv : a) {
        Monome m = Monome(kv.second, kv.first);
        add(ret, m);
    }
    for (auto &kv : b) {
        Monome m = Monome(kv.second, kv.first);
        add(ret, m);
    }
    return ret;
}

void operator+=(Polynome&me,double&d){
    multiset<int> key;
    me[key]+=d;
}
void operator+=(Polynome&p, Monome&s){
    p[s.second]+=s.first;
}
