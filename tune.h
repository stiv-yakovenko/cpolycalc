#ifndef CPOLYCALC_TUNE_H
#define CPOLYCALC_TUNE_H

#include <vector>
#include "interval.h"
#include "Reduced.h"

enum TuneResult {
    NO_ROOT, NO_CHANGE, CHANGE, SPLIT
};

TuneResult tuneQuadratic(vector<interval> &iv, int x, int y, double k/*, vector<double> &sol*/);
TuneResult tuneLinear(Linear &eq, vector<interval> &srcInt);
TuneResult tuneLinearSubst(Subst&subst, vector<interval> &srcInt,
	Reduced&reduced, set<int>&varsToSplit/*,Point&sol*/);


#endif //CPOLYCALC_TUNE_H
