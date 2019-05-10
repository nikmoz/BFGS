#pragma once
struct Matrix{
	double H11 = 0;
	double H12 = 0;
	double H21 = 0;
	double H22 = 0;
};
Matrix SetMatrix(double a11, double a12, double a21, double a22);
Matrix MultiplyMatrixs(Matrix Matrix1, Matrix Matrix2);
Matrix MultiplyMatrix(Matrix Matrix, double coeff);
Matrix AddMatrixs(Matrix Matrix1, Matrix Matrix2);
Matrix SubstMatrixs(Matrix Matrix1, Matrix Matrix2);
Matrix InvertMatrix(Matrix Matrix);