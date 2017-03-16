/*double findNonzeroRadius(Reduced&r, Point&p) {
    Polynome t = decompose(r.grandPoly,p);
    cout << "p="<<p<<", tailor="<<t<<endl;
    double c[5];
    double nrm[4];
    for (auto&kv : t) {
        if (kv.first.size()) {
            c[kv.first.size()] += fabs(kv.second);
        } else {
            c[0] = -fabs(kv.second);
        }
    }
    nrm[3] = c[0] / c[4];
    nrm[2] = c[1] / c[4];
    nrm[1] = c[2] / c[4];
    nrm[0] = c[3] / c[4];
    double sol[5];
    cout << c[0] <<" "<< c[1]<<"*r " <<c[2] << "*r*r " << c[3] << "*r*r*r" << c[4] << "*r*r*r*r"<<endl;
    int s = SolveP4(sol,c[3] / c[4],c[2] / c[4],c[1] / c[4],c[0] / c[4]);
    double bestR = DBL_MAX;
    for (int i = 0; i < s; i++) {
        if (sol[i] >= 0 && sol[i] < bestR) {
            bestR = sol[i];
        }
    }
    cout<<"r="<<bestR<<endl;
    return bestR;
}
*/