//
// Created by steve on 10.11.2015.
//
#include <iostream>
#include "print.h"
#include "interval.h"
double EQ_EPS = 1E-12;

interval addInterval(interval i1, interval i2) {
    return interval(i1.first + i2.first, i1.second + i2.second);
};

interval operator+(interval i1, interval i2) {
    return addInterval(i1, i2);
}

interval subInterval(interval i1, interval i2) {
    return interval(i1.first - i2.first, i1.second - i2.second);
};

interval operator-(interval i1, interval i2) {
    return subInterval(i1, i2);
}

interval addConst(interval i, double c) {
    return interval(i.first + c, i.second + c);
};

bool eq(double v1, double v2) {
    return fabs(v1 - v2) <= max(fabs(v1), fabs(v2)) * EQ_EPS;
};

double limitedDiv(double a, double b, interval limitInterval) {
    if (a >= b * limitInterval.second) return limitInterval.second;
    if (a <= b * limitInterval.first) return limitInterval.first;
    return a / b;
};

bool equal(interval i1, interval i2) {
    return eq(i1.first, i2.first) && eq(i1.second, i2.second);
};

interval mult(interval i1, interval i2) {
    double v1 = i1.first * i2.first;
    double v2 = i1.second * i2.second;
    double v3 = i1.first * i2.second;
    double v4 = i1.second * i2.first;
    return interval(min(min(v1, v2), min(v3, v4)), max(max(v1, v2), max(v3, v4)));
};

interval square(interval i) {
    auto v1 = i.first * i.first;
    auto v2 = i.second * i.second;
    if (i.first * i.second >= 0)
        return interval(min(v1, v2), max(v1, v2));
    return interval(0, max(v1, v2));
};

interval sqrtInterval(interval i) {
    return interval(sqrt(max(i.first, 0.)), sqrt(max(i.second, 0.)));
};

interval multiply(interval i, double c) {
    if (c >= 0) return interval(c * i.first, c * i.second);
    return interval(c * i.second, c * i.first);
};

interval narrow(interval &old, interval &n) {
	interval r = old;
	if (n.first > r.first + EQ_EPS) r.first = n.first;
	if (n.second < r.second - EQ_EPS) r.second = n.second;
	return r;
}

double epsMax(double a, double b) {
    if (fabs(a - b) < EQ_EPS) return min(a, b);
    return max(a, b);
};

double epsMin(double a, double b) {
    if (fabs(a - b) < EQ_EPS) return max(a, b);
    return min(a, b);
};
interval EMPTY(1, -1);

bool isEmpty(interval i) {
    return i.second < i.first;
}

interval join(interval i1, interval i2) {
    if (isEmpty(i1)) return i2;
    if (isEmpty(i2)) return i1;
    return interval(min(i1.first, i2.first), max(i1.second, i2.second));
};

interval intersect(interval i1, interval i2) {
    auto ret = interval(epsMax(i1.first, i2.first), epsMin(i1.second, i2.second));
    if (ret.first > ret.second) return EMPTY;
    return ret;
};

interval reduceInterval(interval &i, interval &src) {
    //    cout<<"i="<<i<<endl;
    //  cout<<"src="<<src<<endl;
    //cout<<(i.first < src.first)<<endl;
    auto x1 = i.first < src.first ? src.first : (i.first - src.first > EQ_EPS ? i.first : src.first);
    auto x2 = i.second > src.second ? src.second : (src.second - i.second > EQ_EPS ? i.second : src.second);
    //cout <<"("<<x1<<","<<x2<<")"<<endl;
    return interval(x1, x2);
}

interval sort(interval src) {
    return src.first <= src.second ? src : interval(src.second, src.first);
}
// returns interval of two roots, EMPTY if no roots, .first=.second if 1 root

interval solveSquare(double a, double b, double c) {
    auto d = b * b - 4 * a * c;
    if (d < 0) return EMPTY;
    if (d == 0) {
        auto x = -.5 * b / a;
        return interval(x, x);
    }
    interval ret;
    // d > 0
    auto sqD = sqrt(d);
    auto x1 = 0., x2 = 0.;
    if (fabs(a) > fabs(b)) {
        x1 = (-b + sqD) / (a * 2);
        x2 = (-b - sqD) / (a * 2);
    } else { // a < b
        x1 = b > 0 ? (2 * c / (-b - sqD)) : (2 * c / (-b + sqD));
        x2 = -b / a - x1;
    }
    return sort(interval(x1, x2));
};

double evalSq(double a, double b, double x) {
    return a * x * x + b * x;
};
// c is interval
// y = at*t+b*b/4/a, t = x+b/2a

interval quadInt(double a, double b, interval xInt) {
    auto y1 = evalSq(a, b, xInt.first);
    auto y2 = evalSq(a, b, xInt.second);
    if (2 * fabs(a) * xInt.first < -b * sgn(a) &&
            -b * sgn(a) < 2 * fabs(a) * xInt.second) { // inside, a!=0
        if (a > 0) {
            return interval(-.25 * b * b / a, max(y1, y2));
        } else {
            return interval(min(y1, y2), -.25 * b * b / a);
        }
    } else { //outside
        return interval(min(y1, y2), max(y1, y2));
    }
};
// b x = -c;
// c is an interval!

vector<interval> quadRoot(double a, double b, interval c, interval xInt) {
  //  cout << "a=" << a << ", b=" << b << ", c=" << c << ", x=" << xInt << endl;
    if (a == 0) {
        auto r = multiply(c, -1. / b);
        auto res = reduceInterval(r, xInt);
        auto ret = vector<interval>();
        if (!isEmpty(res))ret.push_back(res);
        return ret;
    }
    auto i0 = solveSquare(a, b, c.first);
    auto i1 = solveSquare(a, b, c.second);
//    cout << "i0=" << i0;
  //  cout << " i1=" << i1 << endl;
    if (isEmpty(i0)) {
        i0 = i1;
        i1 = EMPTY;
    }
    if (isEmpty(i1)) {
        vector<interval> ret;
        auto i0Int = reduceInterval(xInt, i0);
        if (!isEmpty(i0Int))ret.push_back(i0Int);
        return ret;
    }
    // make i0 smaller i1
    if (i1.second - i1.first < i0.second - i0.first) {
        auto z = i1;
        i1 = i0;
        i0 = z;
    }
//    cout << "i0=" << i0;
  //  cout << " i1=" << i1 << endl;
    auto int1 = sort(interval(i1.first, i0.first));
    auto int2 = sort(interval(i0.second, i1.second));
//    cout << "int1=" << int1 << ",int2=" << int2 << endl;
    auto x1 = reduceInterval(int1, xInt);
    auto x2 = reduceInterval(int2, xInt);
    //cout << "x1="<<x1<<",x2="<<x2<<endl;
    auto ret = vector<interval>();
    if (!isEmpty(x1)) ret.push_back(x1);
    if (!isEmpty(x2)) ret.push_back(x2);
    return ret;
};

bool inside(double x, interval &i) {
    return x >= i.first && x <= i.second;
};

double fit(double x, interval i) {
    return min(max(x, i.first), i.second);
};
bool inside(interval &parent, interval &child) {
	return inside(child.first, parent) && inside(child.second, parent);
}
