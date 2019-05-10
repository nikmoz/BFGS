#pragma once
#include <cmath>
#include "Vector.h"
#include <map>
#include <algorithm>
#include <functional>

struct Interval
{
	double LBound = 0;
	double RBound = 0;
};

class ODSearch
{
public:
	std::function<double(double, double)> function;

	Vector S;
	Vector Xk;
	int FuctionCalls = 0;
	std::map<std::pair<double, double>, double> Fuction;

	ODSearch(double x1,double x2,std::function<double (double,double)> function);  

	double FunctionCalc(double x1, double x2);

	Interval SvenSearch(double X, double Step);

	double DichotomySearch(double LBound, double RBound, double accuracy);
	double GoldenRatioSearch(double LBound, double RBound, double accuracy);
	double DSKSearch(double LBound, double RBound, double accuracy);

private:
	double FunctionCalc(double Labmda);
};








