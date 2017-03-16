//
// Created by steve on 12.11.2015.
//

#include <iomanip>
#include <sstream>
#include "evaluate.h"
#include "algebra.h"
#include "print.h"

template<typename T>
std::string to_string(T const &value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}

namespace print {

    string valWithSign(double v) {
        string ret = "";
        ret += v >= 0 ? "+" : "";
        ret += to_string(v);
        return ret;
    };

    string monome(const multiset<int> &vars) {
        string ret = "";
        auto cnt = vars.size();
        for (auto i : vars) {
            cnt--;
            ret += "x";
            ret += to_string(i);
            if (cnt > 0) ret += "*";
        }
        return ret;
    };

    string polyToString(Polynome &poly) {
        string str = "";
        for (auto &kStr : poly) {
            str += valWithSign(kStr.second);
            if (kStr.first.size() > 0)
                str += "*";
            str += monome(kStr.first);
            str += " ";
        }
        return str;
    };

    string monToString(Monome &mon) {
        string str = "";
        str += valWithSign(mon.first);
        if (mon.second.size() > 0)
            str += "*";
        str += monome(mon.second);
        str += " ";
        return str;
    };

    string sys(deque<Linear> &sys, vector<double> &sol) {
        string ret = "";
        for (auto p = 0U; p < sys.size(); p++) {
            Linear lin = sys[p];
            double disc = eval::linear(lin, sol);
            ret += linToString(lin, sol) + "(*" + to_string(disc) + "*)";
            if (p < sys.size() - 1)
                ret += " ,\n";
        }
        return ret;
    };

    string sys(deque<Polynome> &sys, vector<double> &sol) {
        string ret = "";
        for (auto p = 0U; p < sys.size(); p++) {
            Polynome poly = sys[p];
            double disc = eval::poly(poly, sol);
            ret += polyToString(sys[p]) + "(*" + to_string(disc) + "*)";
            if (p < sys.size() - 1)
                ret += " ,\n";
        }
        return ret;
    };

    string sq(int x, int y, double k) {
        return valWithSign(k) + "*x" + to_string(x) + "^2-x" + to_string(y);
    };

    double evalSq(double k, int x, int y, Point&sol) {
        return k * pow(sol[x], 2) - sol[y];
    }

    string linToString(Linear&lin, Point&sol) {
        string ret;
        for (auto &k : lin) {
            auto i = k.first;
            if (i == -1) {
                continue;
            }
            ret += valWithSign(lin[i]) + "*x" + to_string(i) + " ";
        }
        ret += string("(*");
        ret += to_string(eval::linear(lin, sol));
        ret += string("*)");
        return ret;
    }

    string linSysToString(deque<Linear>&sys, Point&sol) {
        string ret = "";
        for (auto l = 0U; l < sys.size(); l++) {
            auto lin = sys[l];
            ret += valWithSign(lin[-1]);
            ret += " ";
            ret += linToString(lin, sol);
            if (l != sys.size() - 1) ret += ", \n";
        }
        return ret;
    }

    string reduced(Reduced&red, Point& sol) {
        string ret = "{";
        for (auto &square : red.squares) {
            double k = square.second;
            int y = square.first.first;
            int x = square.first.second;
            double disc = evalSq(k, x, y, sol);
            ret = ret + sq(x, y, k);
            ret += ",(*" + to_string(disc) + "*) \n";
        }

        ret += linSysToString(red.linears, sol);

        ret += string("}");
        return ret;
    };

    template <typename T>
    std::string to_string_prec(const T a_value) {
        std::ostringstream out;
        out << std::setprecision(20) << a_value;
        return out.str();
    }

    string point(Point&sol) {
        string vars = "{";
        for (auto i = 0U; i < sol.size(); i++) {
            vars += "x";
            vars += to_string(i);
            vars += "->";
            vars += to_string_prec(sol[i]);
            if (i < sol.size() - 1) vars += ",";
        }
        vars += "}";
        return vars;
    }

    string vect(Point&sol) {
        string vars = "{";
        for (auto i = 0U; i < sol.size(); i++) {
            vars += to_string_prec(sol[i]);
            if (i < sol.size() - 1) vars += ",";
        }
        vars += "}";
        return vars;
    }

    string normal(Linear&linear, Point&sol) {
        string vars = "{";
        for (auto i = 0U; i < sol.size(); i++) {
            vars += to_string_prec(linear[i]);
            if (i < sol.size() - 1) vars += ",";
        }
        vars += "}";
        return vars;
    }
};

ostream& operator<<(ostream& stream, interval& iv) {
    stream << "(" << iv.first << "," << iv.second << ")";
    return stream;
}

ostream& operator<<(ostream& stream, Box& iv) {
    stream << "(";
    for (interval &i : iv) {
        stream << i << ", ";
    }
    stream << ")";

    return stream;
}

ostream& operator<<(ostream& stream, Polynome& iv) {
    stream << print::polyToString(iv);
    return stream;
}

ostream& operator<<(ostream& stream, Monome& mon) {
    stream << print::monToString(mon);
    return stream;
}

ostream& operator<<(ostream& stream, Point& iv) {
    stream << "(";
    for (auto &i : iv) {
        stream << i << ", ";
    }
    stream << ")";
    return stream;
}

ostream& operator<<(ostream& stream, vector<Point>& points) {
    for (Point &p : points) {
        stream << p << endl;
    }
    return stream;
}

ostream& operator<<(ostream& stream, const Point& iv) {
    stream << "(";
    for (auto &i : iv) {
        stream << i << ", ";
    }
    stream << ")";
    return stream;
}

ostream& operator<<(ostream& stream, Subst& subst) {
    for (pair<int, double>k : subst.linears) {
        int x = k.first;
        if (x < 0) {
            stream << print::valWithSign(k.second) << " ";
        } else {
            stream << print::valWithSign(k.second) << "*x" << k.first;
        }
    }
    for (pair<int, double>k : subst.squares) {
        stream << print::valWithSign(k.second) << "*x" << k.first << "^2";
    }
    return stream;
}

ostream& operator<<(ostream& stream, multiset<int>& mon) {
    stream << print::monome(mon);
    return stream;
}

ostream& operator<<(ostream& stream, const multiset<int>& mon) {
    stream << print::monome(mon);
    return stream;
}

ostream& operator<<(ostream& stream, Polysys& iv) {
    for (Polynome &p : iv) {
        stream << p << ",";
    }
    stream << ";";
    return stream;
}

ostream& operator<<(ostream&stream, Linear&lin) {
    string ret;
    for (auto &k : lin) {
        auto i = k.first;
        if (i == -1) {
            ret += print::valWithSign(lin[i]) + " ";
            continue;
        }
        ret += print::valWithSign(lin[i]) + "*x" + to_string(i) + " ";
    }
    stream << ret;
    return stream;
}

ostream& operator<<(ostream& stream, deque<Linear>& iv) {
    for (Linear &p : iv) {
        stream << p << ","<<endl;
    }
    return stream;
}

ostream& operator<<(ostream& stream, map<XY, double>&quads) {
    for (auto &p : quads) {
        auto &xy = p.first;
        double&k = p.second;
        stream << "x" << xy.second;
        stream << "=" << print::valWithSign(k) << "*x" << xy.first << "^2,";
    }
    return stream;
}

ostream& operator<<(ostream& stream, Square& sq) {
    stream << "x" << sq.first.first << "=";
    stream << print::valWithSign(sq.second);
    stream << "*";
    stream << "x" << sq.first.second << "^2, ";
    return stream;
}
