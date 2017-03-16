#include <stdlib.h>
#include <algorithm>
#include <chrono>
#include <iostream>

#include "print.h"
#include "tune.h"
#include "calc.h"
#include "Reduced.h"
#include "wolfram.h"
#include "LPSolve.h"
#include "evaluate.h"
#include "neighborhood.h"

using namespace std;

double calcD(Box&box) {
	double d = DBL_MAX;
	for (auto&kv : box) {
		d = min(d, kv.second - kv.first);
	}
	return d;
}

double sqrt1(double x) {
	return x < 0 ? 0 : sqrt(x);
}
Point sol;
void testIneq(FastLinear ineq, Box&box) {
	double d = eval::linear(ineq, sol);
	bool r = d >= 0;
	if (!r&&inside(box, sol)) {
		cout << "PROBLEM! d=" << d << endl;
		getchar();
	}
}
void genInequations(Reduced&r, Box&box, deque<FastLinear>&ret) {
	//	vector<Linear> neq1;
	for (auto sq : r.squares) {
		FastLinear eq;
		int y = sq.first.first;
		int x = sq.first.second;
		double k = sq.second;
		eq[y] = k;
		testIneq(eq, box);
		ret.push_back(eq);
	}
	//vector<Linear> neq2;
/*	for (auto sq : r.squares) {
		int y = sq.first.first;
		int x = sq.first.second;
		double k = sq.second;
		if (box[x].first*box[x].second >= 0) {
			FastLinear eq;
			eq[y] = 1;
			eq[-1] = -k*min(pow(box[x].first, 2), pow(box[x].second, 2));
			testIneq(eq,box);
			ret.push_back(eq);
			eq = FastLinear();
			eq[y] = -1;
			eq[-1] = k*max(pow(box[x].first, 2), pow(box[x].second, 2));
			testIneq(eq,box);
			ret.push_back(eq);
		} else {
			FastLinear eq;
			eq[y] = -1;
			eq[-1] = k*max(pow(box[x].first, 2), pow(box[x].second, 2));
			testIneq(eq,box);
			ret.push_back(eq);
		}
	}*/

	//vector<Linear> neq3;
/*	for (auto sq : r.squares) {
		int y = sq.first.first;
		int x = sq.first.second;
		double k = sq.second;
		if (box[x].first*box[x].second > 0) {
			FastLinear eq;
			eq[x] = sgn(box[x].first);
			eq[-1] = -sqrt1(min(box[y].first, box[y].second) / k);
			testIneq(eq,box);
			ret.push_back(eq);
			eq = FastLinear();
			eq[x] = -sgn(box[x].first);
			eq[-1] = sqrt1(max(box[y].first, box[y].second) / k);
			testIneq(eq,box);
			ret.push_back(eq);
		} else if (box[x].first < 0 && box[x].second < 0) {
			FastLinear eq;
			eq[x] = -1;
			eq[-1] = -sqrt1(min(box[y].first, box[y].second) / k);
			testIneq(eq,box);
			ret.push_back(eq);
			eq = FastLinear();
			eq[x] = 1;
			eq[-1] = sqrt1(max(box[y].first, box[y].second) / k);
			testIneq(eq,box);
			ret.push_back(eq);
		} else if (box[x].first * box[x].second <= 0) {
			FastLinear eq;
			eq[x] = 1;
			eq[-1] = sqrt1(max(box[y].first, box[y].second) / k);
			testIneq(eq,box);
			ret.push_back(eq);
			eq = FastLinear();
			eq[x] = -1;
			eq[-1] = sqrt1(max(box[y].first, box[y].second) / k);
			testIneq(eq,box);
			ret.push_back(eq);
		}
	}*/

	//vector<Linear> neq4;
	for (auto sq : r.squares) {
		int y = sq.first.first;
		int x = sq.first.second;
		double k = sq.second;
		FastLinear eq;
		eq[y] = 1;
		eq[x] = -2.*k*box[x].second;
		eq[-1] = k*box[x].second* box[x].second;
		testIneq(eq, box);
		ret.push_back(eq);
		eq = FastLinear();
		eq[y] = 1;
		eq[x] = -2.*k*box[x].first;
		eq[-1] = k*box[x].first* box[x].first;
		testIneq(eq, box);
		ret.push_back(eq);
	}
	//vector<Linear> neq5;
	for (auto sq : r.squares) {
		int y = sq.first.first;
		int x = sq.first.second;
		double k = sq.second;
		FastLinear eq;
		eq[y] = -1;
		eq[x] = k*(box[x].second + box[x].first);
		eq[-1] = -k*box[x].second * box[x].first;
		testIneq(eq, box);
		ret.push_back(eq);
		eq = FastLinear();
		eq[y] = 1;
		eq[x] = -k*(box[x].second + box[x].first);
		eq[-1] = .25*k*(box[x].second + box[x].first)*(box[x].second + box[x].first);
		testIneq(eq, box);
		ret.push_back(eq);
	}
}
void solve(Reduced&r, Box box, int depth, double percent) {
	int cnt = 0;
	while (true) {
		cout << "\r";
		cout << floor(percent * 100) << "%, R=  " << getBoxR(box) << ", inside=" << inside(box, sol) << "               ";
		//	cout <<"boxR="<<getBoxR(box)<<endl;
			//cout << "box=" << box << endl;
		cnt++;
		Box prevBox = box;
		for (int z = 0; z < 10; z++) {
			for (auto&pair : r.squares) {
				int var = pair.first.second;
				//int var = floor(rand()*r.numVars/RAND_MAX);
				auto boxV = calcV(box);
				cout << "V=" << boxV << endl;
				auto d = calcD(box);
				if (d < 1E-2) {
					cout << "solution=" << box << endl;
					return;
				}
				for (auto&linear : r.linears) tuneLinear(linear, box);
				for (auto&sq : r.squares) tuneQuadratic(box, sq.first.second, sq.first.first, sq.second);
				deque<FastLinear>ineqs, linears;
				genInequations(r, box, ineqs);
				bool wasinside = inside(box, sol);
				for (auto&linear : r.linears) {
					FastLinear fl;
					for (auto&kv : linear) {
						fl[kv.first] = kv.second;
					}
					linears.push_back(fl);
				}
				auto newInt = solveInterval(var, box, ineqs, linears, r.numVars);
				/*			auto l1 = box[var].second-box[var].first;
							auto l2 = newInt.second-newInt.first;
							cout << l2 / l1 << endl;*/
				if (isEmpty(newInt)) {
					if (wasinside) {
						cout << "Problem" << endl;
						getchar();
					}
					return;
				}
				//cout << "newInt=" << newInt << endl;
				box[var] = intersect(box[var], newInt);
			}
			if (inside(prevBox, sol) != inside(box, sol)) {
				cout << "PROBLEM" << endl;
				getchar();
			}
		}
		if (box == prevBox || cnt > 0) {
			int longest = -1;
			for (auto&pair : r.squares) {
				int x = pair.first.second;
				auto d = box[x].second - box[x].first;
				if (longest < 0 || d >(box[longest].second - box[longest].first)) {// 
					longest = x;
				}
			}
			auto in = box[longest];
			auto c = 0.5*(in.first + in.second);
			box[longest] = interval(in.first, c);
			solve(r, box, depth + 1, percent);
			box[longest] = interval(c, in.second);
			solve(r, box, depth + 1, percent + pow(2, -depth));
			return;
		}
	}
}

void solve2(Reduced&r, Box box, int depth, double percent) {
	int cnt = 0;
	while (true) {
		cout << "\r";
		cout << floor(percent * 100) << "%, R=  " << getBoxR(box) << ", inside=" << inside(box, sol) << "               ";
		//	cout <<"boxR="<<getBoxR(box)<<endl;
		//cout << "box=" << box << endl;
		cnt++;
		Box prevBox = box;

		int longest = -1;
		for (auto&pair : r.squares) {
			int x = pair.first.second;
			auto d = box[x].second - box[x].first;
			if (longest < 0 || d >(box[longest].second - box[longest].first)) {// 
				longest = x;
			}
		}
		int var = longest;
		auto boxV = calcV(box);
		cout << "V=" << boxV << endl;
		auto d = calcD(box);
		if (d < 1E-2) {
			cout << "solution=" << box << endl;
			return;
		}
		for (auto&linear : r.linears) tuneLinear(linear, box);
		for (auto&sq : r.squares) tuneQuadratic(box, sq.first.second, sq.first.first, sq.second);
		deque<FastLinear>ineqs, linears;
		genInequations(r, box, ineqs);
		bool wasinside = inside(box, sol);
		for (auto&linear : r.linears) {
			FastLinear fl;
			for (auto&kv : linear) {
				fl[kv.first] = kv.second;
			}
			linears.push_back(fl);
		}
		auto newInt = solveInterval(var, box, ineqs, linears, r.numVars);
		/*			auto l1 = box[var].second-box[var].first;
		auto l2 = newInt.second-newInt.first;
		cout << l2 / l1 << endl;*/
		if (isEmpty(newInt)) {
			if (wasinside) {
				cout << "Problem" << endl;
				getchar();
			}
			return;
		}
		//cout << "newInt=" << newInt << endl;
		box[var] = intersect(box[var], newInt);
		for (auto&linear : r.linears) tuneLinear(linear, box);
		for (auto&sq : r.squares) tuneQuadratic(box, sq.first.second, sq.first.first, sq.second);

		auto in = box[longest];
		auto c = 0.5*(in.first + in.second);
		box[longest] = interval(in.first, c);
		solve2(r, box, depth + 1, percent);
		box[longest] = interval(c, in.second);
		solve2(r, box, depth + 1, percent + pow(2, -depth));
		return;

	}
}

void main() {
	chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	//	int precision = std::numeric_limits<double>::max_digits10;
	cout.precision(7);
	//	string str = string("List(-Power(x0,2) + x1,-3*Power(x2,2) + x3,-Power(x3,2) + x4,x0 + x1 - x3,x1 + x2 + x3)");
	string str = string("List(-4 + 2*x0 - 2*x0*x4 + 2*x2*x4 + 2*x5 - 2*x0*x5,") +
		"-2 + 2*x1 - 2*x1*x4 + 2*x3*x4 + 4*x5 - 2*x1*x5," +
		"-6 + 2*x2 + 2*x0*x4 - 2*x2*x4,-2 + 2*x3 + 2*x1*x4 - 2*x3*x4," +
		"1 - Power(x0,2) - Power(x1,2) + 2*x0*x2 - Power(x2,2) + 2*x1*x3 -" +
		"Power(x3,2),-4 + 2*x0 - Power(x0,2) + 4*x1 - Power(x1,2))";

	auto cvt = parseWolfram(str);
	Reduced r(cvt);
	Point s0 = { 1.854287547490666901948433658293144957452,		1.480199282222137035802759984080942060162,
		0.93201831765686447855574078078646330060,		1.866747426315165791846532788419724840895,
		2.242275482524551179839123394499501164987,		-2.59127519926380905722981617157019568149 };
	sol = r.updateNewVars(s0);
	cout << "sol=" << sol << endl;
	Box box;
	for (auto i = 0; i < r.numVars; i++) {
		box.push_back(interval(-100, 100));
	}
	cout << "inside=" << inside(box, sol) << endl;
	int var = 0;
	int cnt = 0;
	solve2(r, box, 0, 0);

	chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec" << endl;
	cout << "done!" << endl;
	getchar();
}