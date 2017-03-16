#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace std;
#include "quartic.h"
#include "cubic.h"

vector<double> doSolveQuartic(vector<double>&a,double x0,double x1){
    vector<double> roots;
    return roots;
}
vector<double> solveQuartic(vector<double>a){
    auto s=0.;
    for(auto &e:a){
        s+=fabs(e);
    }
    double R = max(1.,s/fabs(a[a.size()-1]));
    vector<double> q({a[1]/4.,2*a[2]/4.,3*a[3]/4.});
    auto minmax =solveCubic(q[0],q[1],q[2]);
    cout << minmax[0]<<endl;
    cout << minmax[1]<<endl;
    cout << minmax[2]<<endl;
    return doSolveQuartic(a,-R,R);
}
