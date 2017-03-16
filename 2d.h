#ifndef H2D_H
#define H2D_H

#include <vector>

typedef std::pair<double, double> Point2D;
Point2D operator-(Point2D &p1, Point2D&p2); 
Point2D operator+(Point2D &p1, Point2D&p2);
Point2D operator/(Point2D &p1, double c);
double operator*(Point2D&p1, Point2D&p2);
Point2D operator*(Point2D&p1, double c);

#endif