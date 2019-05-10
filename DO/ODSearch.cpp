#include "ODSearch.h"

ODSearch::ODSearch(double x1,double x2,std::function<double(double, double)> function) :function(function) {
	Xk = SetVector(x1, x2);
};

double ODSearch::FunctionCalc(double x1, double x2) {
	double F = function(x1,x2);
	for (auto x : Fuction) {
		if ((x.first.first == x1) && (x.first.second == x2)) {
			return x.second;
		}
	}
	std::pair<double, double> X = { x1,x2 };
	Fuction.emplace(X, F);
	FuctionCalls++;
	return F;
}
double ODSearch::FunctionCalc(double Labmda) {
	Vector X = AddVector(Xk, MultiplyVector(Labmda, S));
	return FunctionCalc(X.x1, X.x2);
}


Interval ODSearch::SvenSearch(double X, double Step) {

	int k = 1;
	if (FunctionCalc(X + Step) > FunctionCalc(X - Step)) {
		Step *= -1;
	}
	X += Step;

	while (FunctionCalc(X) > FunctionCalc(X + pow(2, k) * Step)) {
		X += pow(2, k) * Step;
		k++;
	}
	X += pow(2, k - 1) * Step;

	Interval tmp;
	if (FunctionCalc(X) > FunctionCalc(X - pow(2, k - 1) * Step)) {
		tmp.LBound = X - pow(2, k) * Step;
		tmp.RBound = X;
	}
	else
	{
		tmp.LBound = X - pow(2, k - 1) * Step;
		tmp.RBound = X + pow(2, k - 1) * Step;
	}
	return tmp;
}

double ODSearch::DichotomySearch(double LBound, double RBound, double accuracy) {
	double Middle = (LBound + RBound) / 2;

	while ((RBound - LBound) > accuracy) {
		double LMiddle = (LBound + Middle) / 2;
		double RMiddle = (RBound + Middle) / 2;

		if (FunctionCalc(LMiddle) < FunctionCalc(Middle)) {
			RBound = Middle;
			Middle = LMiddle;
		}

		if (FunctionCalc(RMiddle) < FunctionCalc(Middle)) {
			LBound = Middle;
			Middle = RMiddle;
		}
		else {
			LBound = LMiddle;
			RBound = RMiddle;
		}

	}
	return (RBound + LBound) / 2;
}
double ODSearch::GoldenRatioSearch(double LBound, double RBound, double accuracy){

	double x1 = LBound + 0.382 * (RBound - LBound);
	double x2 = LBound + 0.618 * (RBound - LBound);
	double f1 = FunctionCalc(x1);
	double f2 = FunctionCalc(x2);

	do {
		if (f1 <= f2) {
			RBound = x2;
			x2 = x1;
			f2 = f1;
			x1 = LBound + 0.382 * (RBound - LBound);
			f1 = FunctionCalc(x1);
		}
		else {
			LBound = x1;
			x1 = x2;
			f1 = f2;
			x2 = LBound + 0.618 * (RBound - LBound);
			f2 = FunctionCalc(x2);
		}
	} while ((RBound - LBound) > accuracy);

	return (RBound + LBound) / 2;
}
double ODSearch::DSKSearch(double LBound, double RBound, double accuracy) {
	double Middle = (LBound + RBound) / 2;
	double f1 = FunctionCalc(LBound);
	double f2 = FunctionCalc(Middle);
	double f3 = FunctionCalc(RBound);
	double X = Middle + (Middle - LBound) * (f1 - f3) / (2 * (f1 - 2 * f2 + f3));
	double F = FunctionCalc(X);

	while ((abs(Middle - X) > accuracy) && (abs(f2 - F) > accuracy)) {

		RBound = Middle;
		f3 = f2;
		Middle = X;
		f2 = F;

		double a1 = (f2 - f1) / (Middle - LBound);
		double a2 = ((f3 - f1) / (RBound - LBound) - a1) / (RBound - Middle);
		X = (LBound + Middle) / 2 - a1 / (2 * a2);
		F = FunctionCalc(X);
	}
	return X;
}