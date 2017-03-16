//
// Created by steve on 11.11.2015.
//

#include <algorithm>
#include "util.h"

using namespace std;

// even C++11 standard library is a crap...
string trim(string &str) {
    size_t first = str.find_first_not_of(' ');
    if (first == string::npos)
        return "";
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

string join(vector<string> arr, string glue) {
    auto res = string();
    for (auto i = 0U; i < arr.size(); i++) {
        auto part = arr[i];
        res += part;
        if (i != arr.size() - 1) {
            res += glue;
        }
    }
    return res;
}

vector<string> split(const string s, string sep) {
    vector<string> ret;
    long start = 0;
    auto end = s.find(sep);
    while (end != std::string::npos) {
        ret.push_back(s.substr(start, end - start));
        start = end + sep.length();
        end = s.find(sep, start);
    }
    ret.push_back(s.substr(start, end));
    return ret;
}

double minElem(vector<double>&vect) {
	double ret = vect[0];
	for (double &x : vect) {
		ret = min(ret, x);
	}
	return ret;
}

double maxElem(vector<double>&vect) {
	double ret = vect[0];
	for (double &x : vect) {
		ret = max(ret, x);
	}
	return ret;
}
