#include "ODSearch.h"
#include "Matrix.h"
#include <iostream>
Vector CalcGradient(Vector X,double Step,ODSearch& SearchFunction) {
	Vector tmp;
	tmp.x1 = (SearchFunction.FunctionCalc(X.x1 + Step, X.x2) - SearchFunction.FunctionCalc(X.x1,X.x2))/(Step);
	tmp.x2 = (SearchFunction.FunctionCalc(X.x1, X.x2 + Step) - SearchFunction.FunctionCalc(X.x1, X.x2)) / (Step);
	return tmp;
}
double FindLambda(ODSearch& SearchFunction,double accuaracy) {
	double Lambda = 0;
	double Step = 0.1000001;

	Interval tmp=SearchFunction.SvenSearch(Lambda,Step);
	Lambda=SearchFunction.GoldenRatioSearch(tmp.LBound,tmp.RBound, 0.000001);

	return Lambda;
}
Matrix CalcHessianApro(Matrix HessianApro,Vector Sk,Vector Yk) {
	Matrix tmp;

	Matrix I = SetMatrix(1, 0, 0, 1);
	double pk = 1/(Yk.x1*Sk.x1+Yk.x2*Sk.x2);

	Matrix G = MultiplyMatrix(SetMatrix(Sk.x1*Yk.x1,Sk.x1*Yk.x2,Sk.x2*Yk.x1,Sk.x2*Yk.x2),pk);
	Matrix Coeff1 = SubstMatrixs(I,G);

	G = MultiplyMatrix(SetMatrix(Yk.x1 * Sk.x1, Yk.x1 * Sk.x2, Yk.x2 * Sk.x1, Yk.x2 * Sk.x2), pk);
	Matrix Coeff2 = SubstMatrixs(I, G);

	Matrix Coeff3= MultiplyMatrix(SetMatrix(Sk.x1 * Sk.x1, Sk.x1 * Sk.x2, Sk.x2 * Sk.x1, Sk.x2 * Sk.x2), pk);

	tmp = AddMatrixs(MultiplyMatrixs(MultiplyMatrixs(Coeff1,HessianApro),Coeff2),Coeff3);
	return tmp;
}
int main() {
	int IterCount = 0;
	double accuracy = 0.001;
	double Step = 0.00001;
	Matrix HessianApr = SetMatrix(1,0,0,1);
	//2 * pow(x1 - 11, 2) - x1 * x2 + 5 * pow(x2, 2); Nice test function
	//pow(1 - x1, 2) + 100 * pow(x2 - pow(x1, 2), 2); Rozenbrok function
	//x1*x1-x1*x2+x2*x2+9*x1-6*x2+20; Habr function
	ODSearch SearchFunction = ODSearch(11, 11, [](double x1, double x2) {return pow(1 - x1, 2) + 100 * pow(x2 - pow(x1, 2), 2); });

	Vector Gradient = CalcGradient(SearchFunction.Xk, Step, SearchFunction);

	while (true){
		IterCount++;
		SearchFunction.S.x1 = -HessianApr.H11 * Gradient.x1 - HessianApr.H12 * Gradient.x2;
		SearchFunction.S.x2 = -HessianApr.H21 * Gradient.x1 - HessianApr.H22 * Gradient.x2;

		double Lambda = FindLambda(SearchFunction,accuracy);
		//if (Lambda < accuracy) {
		//	HessianApr = SetMatrix(1, 0, 0, 1);
		//	SearchFunction.S.x1 = -HessianApr.H11 * Gradient.x1 - HessianApr.H12 * Gradient.x2;
		//	SearchFunction.S.x2 = -HessianApr.H21 * Gradient.x1 - HessianApr.H22 * Gradient.x2;
		//	Lambda = FindLambda(SearchFunction, accuracy);
		//}

		Vector PrevX = SearchFunction.Xk;
		SearchFunction.Xk = AddVector(SearchFunction.Xk,MultiplyVector(Lambda, SearchFunction.S));
		Vector Sk = SubstVector(SearchFunction.Xk, PrevX);

		Vector PrevGradient = Gradient;
		Gradient = CalcGradient(SearchFunction.Xk, Step, SearchFunction);
		Vector Yk = SubstVector(Gradient, PrevGradient);

		HessianApr = CalcHessianApro(HessianApr, Sk, Yk);

		double CurrentFunction = SearchFunction.FunctionCalc(SearchFunction.Xk.x1,SearchFunction.Xk.x2);
		double PrevFunction = SearchFunction.FunctionCalc(PrevX.x1,PrevX.x2);
		if (((GetVectorNorm(SubstVector(SearchFunction.Xk,PrevX))/ GetVectorNorm(SearchFunction.Xk)) < accuracy)&&(abs(CurrentFunction- PrevFunction)< accuracy)) {
			break;
		}
	};
	std::cout << IterCount << std::endl;
	std::cout << SearchFunction.Xk.x1<<";"<<SearchFunction.Xk.x2 << std::endl;
	std::cout << SearchFunction.FuctionCalls << std::endl;
	std::cout << SearchFunction.FunctionCalc(SearchFunction.Xk.x1, SearchFunction.Xk.x2) << std::endl;
	system("pause");
	return 0;
}