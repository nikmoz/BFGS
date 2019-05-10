#include "Matrix.h"
Matrix SetMatrix(double a11, double a12, double a21, double a22) {
	Matrix tmp;
	tmp.H11 = a11;
	tmp.H12 = a12;
	tmp.H21 = a21;
	tmp.H22 = a22;
	return tmp;
};
Matrix MultiplyMatrixs(Matrix Matrix1, Matrix Matrix2) {
	Matrix tmp;
	tmp.H11 = Matrix1.H11 * Matrix2.H11 + Matrix1.H12 * Matrix2.H21;
	tmp.H12 = Matrix1.H11 * Matrix2.H12 + Matrix1.H12 * Matrix2.H22;
	tmp.H21 = Matrix1.H21 * Matrix2.H11 + Matrix1.H22 * Matrix2.H21;
	tmp.H22 = Matrix1.H21 * Matrix2.H12 + Matrix1.H22 * Matrix2.H22;
	return tmp;
}
Matrix MultiplyMatrix(Matrix Matrix, double coeff) {
	Matrix.H11 = Matrix.H11 * coeff;
	Matrix.H12 = Matrix.H12 * coeff;
	Matrix.H21 = Matrix.H21 * coeff;
	Matrix.H22 = Matrix.H22 * coeff;
	return Matrix;
};
Matrix AddMatrixs(Matrix Matrix1, Matrix Matrix2) {
	Matrix tmp;
	tmp.H11 = Matrix1.H11 + Matrix2.H11;
	tmp.H12 = Matrix1.H12 + Matrix2.H12;
	tmp.H21 = Matrix1.H21 + Matrix2.H21;
	tmp.H22 = Matrix1.H22 + Matrix2.H22;
	return tmp;
}
Matrix SubstMatrixs(Matrix Matrix1, Matrix Matrix2) {
	Matrix tmp;
	tmp.H11 = Matrix1.H11 - Matrix2.H11;
	tmp.H12 = Matrix1.H12 - Matrix2.H12;
	tmp.H21 = Matrix1.H21 - Matrix2.H21;
	tmp.H22 = Matrix1.H22 - Matrix2.H22;
	return tmp;
}
Matrix InvertMatrix(Matrix Matrix) {
	double tmp = Matrix.H11;
	double Determinant = Matrix.H11 * Matrix.H22 - Matrix.H12 * Matrix.H21;
	Matrix.H11 = Matrix.H22 / Determinant;
	Matrix.H22 = tmp / Determinant;
	Matrix.H12 *= -1 / Determinant;
	Matrix.H21 *= -1 / Determinant;
	return Matrix;
}