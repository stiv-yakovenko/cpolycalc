#include "2d.h"

Point2D operator-(Point2D &p1, Point2D&p2) {
	return Point2D(p1.first - p2.first, p1.second - p2.second);
}

Point2D operator+(Point2D &p1, Point2D&p2) {
	return Point2D(p1.first + p2.first, p1.second + p2.second);
}

double operator*(Point2D&p1, Point2D&p2) {
	return p1.first*p2.first + p1.second*p2.second;
}

Point2D operator*(Point2D&p1, double c) {
	return Point2D(p1.first*c, p1.second*c);
}

Point2D operator/(Point2D &p1, double c) {
	return Point2D(p1.first/c, p1.second/c);
}
