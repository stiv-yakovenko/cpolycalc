#include <map>
#include <string>
#include <iostream>
#include <vector>
#include "LPSolve.h"
#include "print.h"
#include "algebra.h"
#include "interval.h"
#include "lpsolve/lp_lib.h"
using namespace std;
void printLPSolve(map<int, double>&ineq,int eq) {
	cout << ineq << (eq == EQ ? "==0;" : ">=0;") << endl;
}
double solveIneq(int varIdx, Box box, deque<FastLinear>ineqs, deque<FastLinear>&eqs, int numVars, bool maximize) {
	lprec*lp = make_lp(0, numVars);
	int *colno = NULL, j, ret = 0;
	REAL *row = NULL;
	for (int i = 0; i < numVars; i++) {
		string colName = string("x") + to_string(i);
		set_col_name(lp, i + 1, (char*)colName.c_str());
	}
	colno = (int *)malloc(numVars * sizeof(*colno));
	row = (REAL *)malloc(numVars * sizeof(*row));
	set_add_rowmode(lp, TRUE);
	auto addEq = [&](FastLinear&ineq, int constr) {
		//printLPSolve(ineq, constr);		
		double freeCoeff = 0;
		int j = 0;
		for (auto&kv : ineq) {
			if (kv.first == -1) {
				freeCoeff = kv.second;
				continue;
			}
			row[j] = kv.second;
			colno[j] = kv.first + 1;
			j++;
		}
		add_constraintex(lp, j, row, colno, constr, -freeCoeff);
	};
	for (auto&ineq : ineqs) addEq(ineq, GE);
	for (auto&eq : eqs) addEq(eq, EQ);
	j = 1;
	for (auto&interval : box) {
		//cout << interval.first << "<x" << (j-1) << "<" << interval.second << ", " << endl;
		set_bounds(lp, j, interval.first, interval.second);
		j++;
	}
	set_add_rowmode(lp, FALSE);
	double *objRow = new double[numVars + 1];
	for (j = 0; j < numVars + 1; j++) objRow[j] = 0;
	objRow[varIdx + 1] = 1;
	set_obj_fn(lp, objRow);
	delete objRow;
	if (maximize) set_maxim(lp); else set_minim(lp);
	set_verbose(lp, CRITICAL);
	//
	ret = solve(lp);
	if (ret == 25) {
		return maximize?box[varIdx].second: box[varIdx].first;
	}
	if (ret) {
		//cout << "NAN ret=" << ret<<endl;
	//	write_lp(lp, "c:\\\\a.txt");
		return NAN;
	}
	get_variables(lp, row);
	//	for (j = 0; j < numVars; j++)
		//	printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
	free(row);
	free(colno);
	auto r = get_objective(lp);
	delete_lp(lp);
	return r;
}
interval solveInterval(int varIdx, Box &box, deque<FastLinear>&ineqs, deque<FastLinear>&eqs, int numVars) {
	auto l = solveIneq(varIdx, box, ineqs, eqs, numVars, false);
	auto r = solveIneq(varIdx, box, ineqs, eqs, numVars, true);
	interval EMPTY(1, -1);
	if (isnan(l) || isnan(r)) {
		return EMPTY;
	}
	return interval(l, r);
}
/*
void main00() {
	deque<map<int, double>>inEqs, eqs;
	map<int, double> eq;
	eq[0] = -1;
	eq[1] = 1;
	inEqs.push_back(eq);
	eq = map<int, double>();
	eq[1] = -1;
	eq[-1] = -1;
	inEqs.push_back(eq);
	eq = map<int, double>();
	eq[0] = 2;
	eq[1] = -1;
	eq[-1] = +4;
	inEqs.push_back(eq);

	eq = map<int, double>();
	eq[0] = 1;
	eq[1] = 1;
	eq[-1] = 6;
	eqs.push_back(eq);

	Box box;
	box.push_back(interval(-5, 5));
	box.push_back(interval(-5, 5));
	auto r = solveInterval(0, box, inEqs, eqs, 2);
	cout << r << endl;
	getchar();
}

*/