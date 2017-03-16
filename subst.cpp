#include <float.h>
#include "Reduced.h"
#include "subst.h"
void doExclude(deque<Linear> &linears, int var) {
	int maxEqIdx = -1;
	double maxCoeff = 0;
	for (unsigned l = 0; l < linears.size(); l++) {
		if (linears[l].count(var) && fabs(linears[l][var]) >= maxCoeff) {
			maxCoeff = fabs(linears[l][var]);
			maxEqIdx = l;
		}
	}
	Linear repl = linears[maxEqIdx];
	linears.erase(linears.begin() + maxEqIdx);
	for (auto &lin : linears) {
		for (auto &kv : repl) {
			if (kv.first == var)
				lin.erase(var);
			else
				lin[kv.first] -= repl[kv.first] / repl[var];
		}
	}
}
/*Subst mkSubst(Reduced&r, int mainVar) {
	set<int> exclude;
	for (int i = 0; i < r.numVars; i++) {
		exclude.insert(i);
	}
	for (auto&kv : r.squares) {
		if (kv.first.first == mainVar || kv.first.second == mainVar) {
			exclude.erase(kv.first.first);
			exclude.erase(kv.first.second);
		}
	}
	deque<Linear> lins = r.linears;
	while (lins.size() > 1) {
		set<int>::iterator i = exclude.begin();
		doExclude(lins, *i);
	}
}

*/