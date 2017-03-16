#include <cstdlib>
#include <limits>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "algebra.h"
#include "parser.h"
#include "print.h"
#include "square.h"
#include "gradient.h"
#include "evaluate.h"
#include "normalize.h"
#include "neighborhood.h"
#include "wolfram.h"
#include "iterative.h"
#include "cubic.h"
#include "quartic.h"
#include "poly34.h"
#include "tune.h"


bool inside(Box&parent, Box&child) {
	for (unsigned i = 0; i < parent.size(); i++) {
		if (!inside(parent[i], child[i])) return false;
	}
	return true;
}

using namespace std;
void tuneSimple(Reduced&r, Box&box) {
	while (true) {
		Box prevBox = box;
		for (auto&linear : r.linears) tuneLinear(linear, box);
		for (auto&sq : r.squares) tuneQuadratic(box, sq.first.second, sq.first.first, sq.second);
		if (!inside(prevBox, box)) {
			cout << "problem with tuneSimple" << endl;
			getchar();
		}
		if (prevBox == box) return;
	}
}
string boxToString(Box&b) {
	stringstream r;
	for (auto&kv : b) {
		uint64_t *b1 = (uint64_t*)(&kv.first);
		uint64_t *b2 = (uint64_t*)(&kv.second);
		r << *b1;
		r << " ";
		r << *b2;
		r << " ";
	}
	return r.str();
}
double calcVol(Box&box) {
	double r = 1;
	for (auto i : box) {
		r *= i.second - i.first;
	}
	return r;
}
void boxIterate(Reduced&r, Box box, vector<Point>&allRoots, int depth) {
	/*    cout.precision(17);
		Point sl = {1.7634461203153509072715803028927950560639240811099,
			1.3541284792039218991442153033803225083293561086555,
			2.7247999726431325494720518840050925341568038468,
			1.0788127139187973220560762548608949507249579816996,
			-0.28626298910698669344619665590600085824200004863,
			-0.6703209217045077554834247470797184959849024162400};
		Point s = r.updateNewVars(sl);*/
	int i = 0;
	set<Box> prevBoxes;
	while (true) {
		i++;
		if (prevBoxes.count(box)) {
			std::ofstream out("problem1.txt");
			out << boxToString(box);
			out.close();
			cout << "PROBLEM!" << endl;
			getchar();
		}
		prevBoxes.insert(box);
		Box prevBox = box;
		cout << "========================" << endl;
		cout << "it=" << i << " depth=" << depth << ", rootCount=" << allRoots.size() << " V=" << calcVol(box);
		cout << " box=" << box << endl;
		Point c = center(box);
		auto R = getBoxR(box);
		auto R1 = findNonzeroRadius(r, c);
		cout << "R_nz=" << R1 << ", R=" << R << endl;
		if (R < 1E-6) {
			cout << "SMALL" << endl;
			allRoots.push_back(c);
			return;
		}
		tuneSimple(r, box);
		for (int longest = 0; longest < r.numVars; longest++) {
			vector <interval> roots = tuneCombination(c, r, box, longest);
			vector <interval> goodRoots;
			for (auto&root : roots) {
				auto inter = narrow(box[longest], root);
				if (isEmpty(inter)) continue;
				goodRoots.push_back(inter);
			}
			//cout << "goodRoots = " << goodRoots << endl;
			switch (goodRoots.size()) {
			case 0:
				return;
			case 1:
				if (goodRoots[0] == box[longest]) continue;
				if (!inside(box[longest], goodRoots[0])) {
					cout << "another problem" << endl;
					getchar();
				}
				box[longest] = goodRoots[0];
				continue;
			case 2: {
				if (goodRoots[0] == goodRoots[1]) {
					cout << goodRoots[0] << endl;
					std::ofstream out("problem.txt");
					out << boxToString(box);
					out.close();
					cout << "PROBLEM!" << endl;
					getchar();
					exit(0);
				}
				box[longest] = goodRoots[0];
				boxIterate(r, box, allRoots, depth + 1);
				box[longest] = goodRoots[1];
				boxIterate(r, box, allRoots, depth + 1);
			}
					return;
			}
		}
		if (box == prevBox) {
			int longest = -1, i = 0;
			for (auto&pair : box) {
				auto d = pair.second - pair.first;
				if (longest < 0 || d >(box[longest].second - box[longest].first)) {
					longest = i;
				}
				i++;
			}
			auto in = box[longest];
			box[longest] = interval(in.first, .5*(in.first + in.second));
			boxIterate(r, box, allRoots, depth + 1);
			box[longest] = interval(.5*(in.first + in.second), in.second);
			boxIterate(r, box, allRoots, depth + 1);
			return;
		}
	}
}

int main00() {
	Reduced r("x1-1+2*x0,x1-x0*x0");
	cout << r.linears << endl;
	cout << r.squares << endl;
	Point center = Point({ 1, 0 });
	Box b = vector<interval>({ interval(0, 2), interval(-2, 2) });
	vector<Point>roots;
	boxIterate(r, b, roots, 0);
	cout << "roots=" << roots;
	return 0;
}

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

Box stringToBox(string str) {
	Box b;
	vector<string> parts = split(str, ' ');
	for (auto i = 0U; i < parts.size(); i += 2) {
		std::stringstream S(parts[i]);
		double x;
		S >> *(uint64_t*)&x;
		std::stringstream S1(parts[i + 1]);
		double y;
		S1 >> *(uint64_t*)&y;
		b.push_back(interval(x, y));
	}
	return b;
}



int main1(int argc, char** argv) {
	int precision = std::numeric_limits<double>::max_digits10;
	cout.precision(precision);
	//	cout << "str=" << bx << endl;
		//string str = string("List(x0*x0-1,x0-x1*x1)");
	string str = string("List(-4 + 2*x0 - 2*x0*x4 + 2*x2*x4 + 2*x5 - 2*x0*x5,") +
		"-2 + 2*x1 - 2*x1*x4 + 2*x3*x4 + 4*x5 - 2*x1*x5," +
		"-6 + 2*x2 + 2*x0*x4 - 2*x2*x4,-2 + 2*x3 + 2*x1*x4 - 2*x3*x4," +
		"1 - Power(x0,2) - Power(x1,2) + 2*x0*x2 - Power(x2,2) + 2*x1*x3 -" +
		"Power(x3,2),-4 + 2*x0 - Power(x0,2) + 4*x1 - Power(x1,2))";
	auto cvt = parseWolfram(str);
	Reduced r(cvt);
	std::ofstream out("sq.txt");
	out << r.toWolfram() << endl;
	out.close();
	Point sol0 = { 1.854287547490666901948433658293144957452,
		1.480199282222137035802759984080942060162,
		0.93201831765686447855574078078646330060,
		1.866747426315165791846532788419724840895,
		2.242275482524551179839123394499501164987,
		-2.59127519926380905722981617157019568149 };
	Point sol = { 1.7634461203153509072715803028927950560639240811099,
		1.3541284792039218991442153033803225083293561086555,
		2.7247999726431325494720518840050925341568038468,
		1.0788127139187973220560762548608949507249579816996,
		-0.28626298910698669344619665590600085824200004863,
		-0.6703209217045077554834247470797184959849024162400 };
	Point s = r.updateNewVars(sol);
		Box b;
		for (auto i = 0; i < r.numVars; i++) {
			b.push_back(interval(s[i] - 100, s[i] + 100));
		}
	//Box b = stringToBox("13722977014762069110 4611686018564350144 4607182418419014134 4612846177504910906 13819895845540795632 4611959207651258666 13719853039567765504 4611686018501885356 13832956733395175091 4612814224450769689 13836425475440147014 13826708157639821832 0 4616189618328682884 4607182418038010876 4618808801383559846 0 4616752568213623074 0 4616189618203753305 0 4618728659317373938 13836186261454931122 4613868655359180957 4616649354760995054 4621152954388365550 4599585756122372817 4619339645785941755 4612976243874292546 4615376124798074613 4619139702539822827 4623643302167193323 13836433618316094119 4613061581461318311 4612059112844270720 4619360891075473593 13832810867892058423 4612887157202328023 4607596200824278800 4618912246984876010 4612926119553500411 4615305238303805636 4619011291780965973 4623514891408336469 13836186261400796316 4614791430153164317 4618544863714005536 4622617482270357283 13830554455654793216 4607182418800017408 0 4607182418800017408 13830554455654793216 4607182418800017408 0 4607182418800017408");
	vector<Point>roots;
	boxIterate(r, b, roots, 0);
	cout << "roots=" << roots;
	getchar();
	return 0;
}

