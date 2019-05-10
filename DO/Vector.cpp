#include "Vector.h"
Vector AddVector(Vector vector1, Vector vector2) {
	Vector tmp;
	tmp.x1 = vector1.x1 + vector2.x1;
	tmp.x2 = vector1.x2 + vector2.x2;
	return tmp;
}
Vector SubstVector(Vector vector1, Vector vector2) {
	Vector tmp;
	tmp.x1 = vector1.x1 - vector2.x1;
	tmp.x2 = vector1.x2 - vector2.x2;
	return tmp;
}
Vector MultiplyVector(double K, Vector vector) {
	Vector tmp;
	tmp.x1 = vector.x1 * K;
	tmp.x2 = vector.x2 * K;
	return tmp;
}
double MultiplyVectors(Vector vector1, Vector vector2) {
	return vector1.x1* vector2.x1 + vector1.x2 * vector2.x2;
}
Vector SetVector(double x1, double x2) {
	Vector tmp;
	tmp.x1 = x1;
	tmp.x2 = x2;
	return tmp;
}
double GetVectorNorm(Vector X) {
	return pow(X.x1 * X.x1 + X.x2 * X.x2, 0.5);
};