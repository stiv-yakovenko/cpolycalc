#ifndef GRADIENT_H
#define GRADIENT_H
	  
#include "algebra.h"

void iterDescent(Reduced&r, Point &p0);
void calcDerivativesFast(Reduced&r);
void descentQuad(Reduced&r, Point &p0);
Point minimizeQForm(Polynome qTailor,Point p);
Polynome quadTailor(Reduced&r,Point &p);
void calcG(Reduced &r);
void calcGradSteps(Reduced&ren);
void minimize(Reduced &reduced, Point&sol);
void calcGDerivatives(Reduced &r);
Polynome tailor(Reduced&r, Point&p);
Polynome dirPoly (Reduced&r,Point&p,Point&v);
void descent(Reduced&r, Point &p0);
#endif