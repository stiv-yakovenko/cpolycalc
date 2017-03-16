#include <math.h>
#include "normalize.h"

double calcK(Linear &n, Linear&m) {
    double ret = 0, m2 = 0.;

    for (auto &key : m) {
        if (key.first == -1) continue;
        m2 += pow(key.second, 2);
        if (n.count(key.first) == 0) continue;
        ret += key.second * n[key.first];
    }
    return ret / m2;
}

void normalize(LinSyst&sys, Point&sol) {
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