//
// Created by steve on 11.11.2015.
//

#ifndef CPOLYCALC_PARSER_H
#define CPOLYCALC_PARSER_H

#include <string>
#include <set>
#include "util.h"
#include "algebra.h"
using namespace std;

Polynome parsePoly(string polyStr);

namespace parser {
    Polysys system(string systemStr);
}

#endif //CPOLYCALC_PARSER_H
