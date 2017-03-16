#include <iostream>
#include "tune.h"
#include "algebra.h"
#include "interval.h"
#include "print.h"

using namespace std;

class Tunable {
public:
	FastLinear *linear = 0;
	pair <XY, double> quad;
};

void smartTune(Reduced&r) {
	vector<Tunable> tunables;
	multimap<double, Tunable> discrepancies;
	for (auto&kv : r.squares) {
		Tunable t;
		t.quad = pair<XY, double>(kv.first, kv.second);
		tunables.push_back(t);
		discrepancies.insert(pair<double,Tunable>(0, t));
	}
	for (auto&linear : r.linears) {
		Tunable t;
		t.linear = new FastLinear();
		FastLinear fl;
		for (auto&kv : linear) {
			fl[kv.first] = kv.second;
		}
		tunables.push_back(t);
		discrepancies.insert(pair<double, Tunable>(0, t));
	}
/*	while (true) {
		d
		
	}*/
	for (auto&tunable : tunables) {
		if (tunable.linear) delete tunable.linear;
	}
}

interval calcSum(int j, vector<interval> &srcInt, Linear &eq) {
	interval sum = interval(0, 0);
	for (auto &key : eq) {
		int i = key.first;
		if (i == j) continue;
		auto iv = (i == -1) ? interval(1, 1) : srcInt[i];
		iv = multiply(iv, eq[i]);
		sum = addInterval(sum, iv);
	}
	return sum;
}

TuneResult tuneLinear(Linear &eq, vector<interval> &srcInt) {
	bool somethingChanged = false;
	for (auto &key : eq) {
		int j = key.first;
		if (j == -1) continue;
		interval sum = calcSum(j, srcInt, eq);
		interval jInt = multiply(sum, -1. / eq[j]);
		interval newJInt = reduceInterval(jInt, srcInt[j]);
		if (isEmpty(newJInt)) return NO_ROOT; // end of story
		somethingChanged = somethingChanged || !equal(srcInt[j], newJInt);
		srcInt[j] = newJInt;
	}
	return somethingChanged ? CHANGE : NO_CHANGE;
};

double defGet(Linear&lin, int idx) {
	if (lin.count(idx) == 0) {
		return 0;
	} else {
		return lin[idx];
	}
}

TuneResult tuneLinearSubst(Subst&subst, vector<interval> &srcInt,
	Reduced&reduced, set<int>&varsToSplit/*, Point&sol*/) {
	set<int> vars;
	for (auto &k : subst.linears) vars.insert(k.first);
	for (auto &k : subst.squares) vars.insert(k.first);
	vars.erase(-1);
	bool somethingChanged = false;
	auto c = defGet(subst.linears, -1);
	for (auto v : vars) {
		interval freePart = interval(c, c);
		for (auto w : vars) {
			if (w == v) continue;
			auto sqInt = quadInt(defGet(subst.squares, w), defGet(subst.linears, w), srcInt[w]);
			freePart = addInterval(freePart, sqInt);
		}
		auto sqRoots = quadRoot(defGet(subst.squares, v), defGet(subst.linears, v), freePart, srcInt[v]);
		if (sqRoots.size() == 0) {
			return NO_ROOT;
		}
		if (sqRoots.size() == 2) {
			sqRoots = { interval(sqRoots[0].first, sqRoots[1].second) };
		}
		auto prev = srcInt[v];
		somethingChanged = somethingChanged || !equal(srcInt[v], prev);
	}
	return somethingChanged ? CHANGE : NO_CHANGE;
};

TuneResult tuneQuadratic(vector<interval> &iv, int x, int y, double k/*, vector<double> &sol*/) {
	bool somethingChanged = false;
	interval r = multiply(square(iv[x]), k);
	interval yInt = reduceInterval(r, iv[y]);
	if (isEmpty(yInt)) return NO_ROOT;
	somethingChanged = somethingChanged || !equal(iv[y], yInt);
	iv[y] = yInt;
	// x = sqrt(y/k)
	interval m = multiply(yInt, 1. / k);
	interval sqrt1 = sqrtInterval(m);
	if (isEmpty(sqrt1)) return NO_ROOT;
	interval sqrt2 = multiply(sqrt1, -1);
	interval x1 = reduceInterval(sqrt1, iv[x]);
	interval x2 = reduceInterval(sqrt2, iv[x]);
	if (isEmpty(x1) && isEmpty(x2)) return NO_ROOT;
	interval xInt = join(x1, x2);
	somethingChanged = somethingChanged || !equal(iv[x], xInt);
	if (isEmpty(xInt)) return NO_ROOT;
	iv[x] = xInt;
	return somethingChanged ? CHANGE : NO_CHANGE;
};

