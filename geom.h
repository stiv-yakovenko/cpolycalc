
#ifndef GEOM_H
#define GEOM_H

Point closestPoint(Point p1,Point p2,Point v1,Point v2);
double dist(Linear&line,Point &p);
double distToQuad(Point &sol, int x, int y, double k);

#endif /* GEOM_H */

