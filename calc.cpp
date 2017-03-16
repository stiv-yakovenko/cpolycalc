//
// Created by steve on 13.11.2015.
//

#include <chrono>
#include <algorithm>
#include <iostream>
#include "calc.h"
#include "algebra.h"
#include "interval.h"
#include "tune.h"
#include "evaluate.h"
#include "print.h"

using namespace std;
using namespace std::chrono;

int doTune(Reduced &reduced, Box &box/*, Point &sol, pair<XY, double>&worstSQ*/) {
	int cnt = 0;
	int worstSqIdx = -1;
	pair<XY, double> ret;
	set<int> varsToSplit;
	double worstSqD;
	while (true) {
		worstSqD = 0;
		cnt++;
		bool somethingChanged = false;
		for (auto key : reduced.squares) {
			//var vars = sqKey.split(',');
			int y = key.first.first;
			int x = key.first.second;
			double k = key.second;
			TuneResult res = tuneQuadratic(box, x, y, k);
			if (res == NO_ROOT) {
/*				if (ins) {
					cout << "Problem" << endl;
				}*/
				return -1;
			}
			double d = box[x].second - box[x].first;
			if (fabs(d) > worstSqD) {
				worstSqIdx = x;
				worstSqD = fabs(d);
				//worstSQ = pair<XY, double>(XY(x, y), k);
			}
			somethingChanged = somethingChanged || (res == CHANGE);
		}
		for (Subst subst : reduced.substs) {
			TuneResult res = tuneLinearSubst(subst, box, reduced, varsToSplit);
			if (res == NO_ROOT) {
//				if (ins) {
//					cout << "Problem" << endl;
//				}
				return -1;
			}
			somethingChanged = somethingChanged || (res == CHANGE);
		}
		if (!somethingChanged) break;
	}

	if (varsToSplit.size() > 0) {
		int x = *varsToSplit.begin();
		auto ky = reduced.xToKY[x];
		//worstSQ = pair<XY, double>(XY(x, ky.first), ky.second);
		return x;
	}

	return worstSqIdx;
};

const double EPS = 1E-5;
auto prevTS = (long)system_clock::now().time_since_epoch().count();
vector<vector<double>> roots;

double dist(vector<double> &sol, vector<interval> &iv) {
	double d = 0;
	for (auto i = 0U; i < sol.size(); i++) {
		d += pow((sol[i] - (iv[i].first + iv[i].second) / 2), 2);
	}
	return sqrt(d);
};

double calcV(vector<interval> &iv) {
	double V = 1.;
	for (auto i = 0U; i < iv.size(); i++) {
		V *= (iv[i].second - iv[i].first);
	}
	return V;
};

/*void addRoot(vector<Point> &roots, Point &root) {
	for (Point&other:roots) {
		if (length(root - other)<EPS*10) return;
	}
	roots.push_back(root);
}*/
// D -- interval diameter, V -- interval volume
bool boxCalc(Reduced &reduced, Box box, int depth, Point &sol, vector<Point>&roots, double progress) {
	//dout << "iterating:" << depth << endl;

	pair<XY, double> worstSq;

	bool before = inside(box, sol);
	int xIdx = doTune(reduced, box/*, sol/*, worstSq*/);
	bool after = inside(box, sol);

	if (before) {
		cout << "inside, depth = " << depth << endl;
	}

	if (before != after) {
		cout << "Problem" << endl;
		fgetc(stdin);
	}

	if (xIdx == -1) {
		return false;
	}

	double maxD = 0;

	for (auto i = 0U; i < box.size(); i++) {
		double d = box[i].second - box[i].first;
		if (d > maxD) {
			maxD = d;
		}
	}

	if (maxD < EPS) {
		vector<double> root;
		for (auto i = 0U; i < box.size(); i++) {
			root.push_back((box[i].second + box[i].first) / 2);
		}
		//addRoot(roots,root);
		//dout << "sol = " << root << endl;
		return true;
	}

	long ts = (long)system_clock::now().time_since_epoch().count();

	//cout << ts << "," << prevTS << endl; 
	if (abs(ts - prevTS) > 20000000) {
		cout << "inside: " << before << endl;
		cout << "depth = " << depth << endl;
		cout << "V = " << calcV(box) << endl;
		cout << "iv = " << box << endl;
		cout << "progress = " << progress << endl;
		cout << "D = " << maxD << endl;
		cout << "discr = " << dist(sol, box) << endl;
		cout << "Root count = " << roots.size();
		cout << endl;
		prevTS = ts;
	}

	double k = worstSq.second;
	interval oldX = box[xIdx];
	double cX;
	if (oldX.first < 0 && oldX.second > 0) {
		cX = 0;
	} else {
		cX = .5*(oldX.first + oldX.second);
	}
	double cY = k * cX*cX;

	int yIdx = worstSq.first.second;
	interval oldY = box[yIdx];
	// left interval
	vector<double> yVals;// (cX > 0) ? {} : {};
	yVals.push_back(k*oldX.first*oldX.first);
	yVals.push_back(cY);
	box[xIdx] = interval(oldX.first, cX);
	if (inside(0, box[xIdx])) {
		yVals.push_back(0);
	}
	auto y1 = interval(minElem(yVals), maxElem(yVals));
	box[yIdx] = y1;
	bool i1 = inside(box, sol);
	boxCalc(reduced, box, depth + 1, sol,roots, progress);

	// second interval
	yVals.clear();
	yVals.push_back(k*oldX.second*oldX.second);
	yVals.push_back(cY);
	box[xIdx] = interval(cX, oldX.second);
	if (inside(0, box[xIdx])) {
		yVals.push_back(0);
	}
	auto y2 = interval(minElem(yVals), maxElem(yVals));
	box[yIdx] = y2;
	bool i2 = inside(box, sol);
	boxCalc(reduced, box, depth + 1, sol,roots, progress + pow(.5,depth+1));

	if (before && !(i1 || i2)) {
		cout << "Problem" << endl;
		fgetc(stdin);
	}

	box[xIdx] = oldX;
	box[yIdx] = oldY;
	return true;
};
