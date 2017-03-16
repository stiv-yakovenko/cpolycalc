#ifndef CPOLYCALC_ALGEBRA_H
#define CPOLYCALC_ALGEBRA_H

#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <deque>

using namespace std;

typedef std::pair<double, std::multiset<int> > Monome;
typedef std::map<std::multiset<int>, double> Polynome;
typedef std::map<int, double> Linear;
typedef std::unordered_map<int, double> FastLinear;
typedef std::deque<Polynome> Polysys;
typedef std::vector<double> Point;
typedef std::vector<double> Vector;
typedef std::deque<Linear> LinSyst;
typedef std::pair<std::pair<int, int>, double> Square;
typedef std::vector<std::pair<double, double> > Box;
typedef std::pair<int, int> XY;

struct Subst {
    map<int, double> linears;
    map<int, double> squares;
};


// x1 -= x2
void subtract(Point &x1, Point &x2);
Point sub(Point&p1, Point&p2);
double length2(Point&p);
double length(Point&p);
double operator*(Point&p1, Point&p2);
Point operator*(Point&p1, double c);
Point operator+(Point&p1, Point&p2);
void operator+=(Point&me,Point&d);
void operator+=(Polynome&me,double&d);
void operator-=(Point&me,Point&d);
Point operator-(Point&p1, Point&p2);
Point operator/(Point&p1, double c);
void operator/=(Polynome&p, double c);
void operator-=(Polynome&p, Polynome&s);
void operator+=(Polynome&p, Polynome&s);
void operator+=(Polynome&p, Monome&s);
void operator*=(Polynome&p, double c);
bool inside(Box&box, Point&p);
void add(Polynome&poly,Monome&mon);
Monome mult(Monome&a,Monome&b);
Polynome mult(Polynome&a,Polynome&b);
Polynome operator*(Polynome&a,Polynome&b);
Polynome operator+(Polynome&a,Polynome&b);
void removeZeroes(Polynome&p);


#endif //CPOLYCALC_ALGEBRA_H
