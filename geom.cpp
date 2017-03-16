#include <iostream>
#include "algebra.h"
#include "math.h"
#include "geom.h"
#include "print.h"
#include "iterative.h"
using namespace std;

double dist(Linear&line, Point &p) {
    double s = 0, q = 0;
    for (auto&kv : line) {
        if (kv.first == -1) {
            s += kv.second;
        } else {
            s += kv.second * p[kv.first];
            q += kv.second * kv.second;
        }
    }
    return fabs(s) / sqrt(q);
}

double distToQuad(Point &sol, int x, int y, double k) {
    Point2D p = projectToQuad(sol,x,y,k);
    return sqrt(pow(p.first-sol[x],2)+pow(p.second-sol[y],2));
}

Point closestPoint(Point p1, Point p2, Point v1, Point v2) {
    Point p = p1 - p2;
    double v1v2 = v1*v2;
    double v1v1 = v1*v1;
    double v2v2 = v2*v2;
    double v2p = v2*p;
    double v1p = v1*p;
    double d = pow(v1v2, 2) - v1v1*v2v2;
    double t1 = -(v1v2 * v2p - v1p * v2v2) / d;
    Point vt1 = v1*t1;
    return p1 + vt1;
}

/*int main(){
    auto p1=Point({1,0,1});
    auto v1=Point({-.1,0,0});
    auto p2=Point({0,0,0});
    auto v2=Point({0,0.1,0.1});
    Point r= closestPoint(p1,p2,v1,v2);
    cout << r<<endl;
}

 */