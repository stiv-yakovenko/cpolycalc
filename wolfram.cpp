//
// Created by steve on 14.11.2015.
//

#include "wolfram.h"
#include <string>
#include <regex>
#include <iostream>

using namespace std;

string genPower(string var, int cnt) {
    string ret = "";
    for (int i = 0; i < cnt; i++) {
        ret += var;
        if (i != cnt - 1) ret += "*";
    }
    return ret;
}


string parseWolfram(string wStr) {
    smatch m;
    string ret = wStr;
    regex rgx("Power\\((.*?),(\\d*?)\\)");
    while (regex_search(ret, m, rgx)) {
        int pwr = atoi(m[2].str().c_str());
        string repl = genPower((string)m[1].str(), pwr);
        long p = m.position(0);
        ret.replace(p, m[0].str().length(), repl);
    }

	ret = regex_replace(ret, regex("List\\("), "");
	ret = regex_replace(ret, regex("\\)"), "");
	ret = regex_replace(ret, regex(","), ";");
	ret = regex_replace(ret, regex(" "), "");

	cout << "After wolfram preparse " << ret << endl;

    return ret;
}