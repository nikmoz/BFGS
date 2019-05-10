#pragma once
#include <cmath>
struct Vector
{
	double x1 = 0;
	double x2 = 0;
};

Vector AddVector(Vector vector1, Vector vector2);
Vector SubstVector(Vector vector1, Vector vector2);
Vector MultiplyVector(double K, Vector vector);
double MultiplyVectors(Vector vector1, Vector vector2);
Vector SetVector(double x1, double x2);
double GetVectorNorm(Vector X);