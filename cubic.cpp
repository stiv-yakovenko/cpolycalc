#include <vector>
#include <math.h>

const double PI = 3.141592653589793238463;
using namespace std;

vector <double> solveCubic(double a, double b, double c) {
    vector<double> ret;
    auto q = (a * a - 3 * b);
    auto r = (2 * a * a * a - 9 * a * b + 27 * c);
    auto Q = q / 9;
    auto R = r / 54;
    auto Q3 = Q * Q * Q;
    auto R2 = R * R;
    auto CR2 = 729 * r * r;
    auto CQ3 = 2916 * q * q * q;
    if (R == 0 && Q == 0) {
        ret.push_back(-a / 3);
        return ret;
    } else if (CR2 == CQ3) {
        auto sqrtQ = sqrt(Q);

        if (R > 0) {
            ret.push_back(-2 * sqrtQ - a / 3);
            ret.push_back(sqrtQ - a / 3);
            return ret;
        } else {
            ret.push_back(-sqrtQ - a / 3);
            ret.push_back(2 * sqrtQ - a / 3);
            return ret;
        }
    } else if (R2 < Q3) {
        auto sgnR = (R >= 0 ? 1 : -1);
        auto ratio = sgnR * sqrt(R2 / Q3);
        auto theta = acos(ratio);
        auto norm = -2 * sqrt(Q);
        ret.push_back(norm * cos(theta / 3) - a / 3);
        ret.push_back(norm * cos((theta + 2.0 * PI) / 3) - a / 3);
        ret.push_back(norm * cos((theta - 2.0 * PI) / 3) - a / 3);
        return ret;
    } else {
        auto sgnR = (R >= 0 ? 1 : -1);
        auto A = -sgnR * pow(fabs(R) + sqrt(R2 - Q3), 1.0 / 3.0);
        auto B = Q / A;
        ret.push_back(A + B - a / 3);
        return ret;
    }
};