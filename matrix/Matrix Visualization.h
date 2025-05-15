//This file contains function declarations of matrix visualization functions.
//Please note that there's no visualizeDeterminant() function. 
//If you need to visualize determinant calculation, call the determinant function with a boolean value true.
#pragma once
#include "matrix.h"
#include <iostream>
namespace MatrixViz {
	//Everything is kept inside the namespace MatrixViz to avoid polluting the global namespace
	//This avoids name clashes with other libraries or user code


	//Displays multiplication of two matrices step by step
	void visualizeMultiplication(const Matrix& mat1, const Matrix& mat2, const Matrix& result);

	//Displays addition of two matrices step by step
	void visualizeAddition(const Matrix& mat1, const Matrix& mat2, const Matrix& result);

	//Displays subtraction of two matrices step by step
	void visualizeSubtraction(const Matrix& mat1, const Matrix& mat2, const Matrix& result);

	//Displays transposition of a matrix step by step
	void visualizeTranspose(const Matrix& mat, const Matrix& result);

	//Displays inverse calculation of a matrix step by step
	void visualizeInverse(const Matrix& mat, const Matrix& inverseMat);

	//Displays LU decomposition of a matrix step by step
	void visualizeLUDecomposition(const Matrix& mat, const Matrix& L, const Matrix& U);

	//Displays QR decomposition of a matrix step by step
	void visualizeQRDecomposition(const Matrix& mat, const Matrix& Q, const Matrix& R);

	//Displays eigenvalue calculation of a matrix step by step
	void visualizeEigenvalues(const Matrix& mat, const std::vector<std::complex<double>>& eigenvalues);

	//Displays Cholesky decomposition of a matrix step by step
	void visualizeCholeskyDecomposition(const Matrix& mat, const Matrix& L);

	//Displays linear equation solution step by step
	void visualizeLinearEquationSolution(const Matrix& A, const std::vector<double>& x, const std::vector<double>& b);
}


//End of the header file
