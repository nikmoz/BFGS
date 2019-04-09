#include "OneDSearch.h"
#include "Vector.h"
#include <iostream>
double FindLabmdaStep(Vector X, Vector direction) {
	double XModule = sqrt(pow(X.x1,2)+pow(X.x2,2));
	double DModule = sqrt(pow(direction.x1, 2) + pow(direction.x2, 2));
	return 0.1*XModule / DModule;
}
int main() {
	Vector X1, X2;
	double accuracy = 0.001;
	std::cout << FunctionCalc(X0.x1, X0.x2) << std::endl;

	Interval tmp=SvenSearch(0,FindLabmdaStep(X0,S0));
	double Lambda = DichotomySearch(tmp.LBound,tmp.RBound, accuracy);
	X0 = AddVector(X0, MultiplyVector(Lambda, S0));
	std::cout << FunctionCalc(X0.x1, X0.x2) << std::endl;

	X1 = X0;

	S0 = SetVector(0, 1);
	tmp = SvenSearch(0, FindLabmdaStep(X0, S0));
	Lambda = GoldenRatioSearch(tmp.LBound, tmp.RBound, accuracy);
	X0 = AddVector(X0, MultiplyVector(Lambda, S0));
	std::cout << FunctionCalc(X0.x1, X0.x2) << std::endl;

	S0 = SetVector(1, 0);
	tmp = SvenSearch(0, FindLabmdaStep(X0, S0));
	Lambda = DSKSearch(tmp.LBound, tmp.RBound, accuracy);
	X0 = AddVector(X0, MultiplyVector(Lambda, S0));
	std::cout << FunctionCalc(X0.x1, X0.x2) << std::endl;

	X2 = X0;

	S0 = AddVector(X2,MultiplyVector(-1,X1));
	tmp = SvenSearch(0, FindLabmdaStep(X0, S0));
	Lambda = DSKSearch(tmp.LBound, tmp.RBound, accuracy);
	X0 = AddVector(X0, MultiplyVector(Lambda, S0));
	std::cout << FunctionCalc(X0.x1, X0.x2) << std::endl;
	return 0;
}