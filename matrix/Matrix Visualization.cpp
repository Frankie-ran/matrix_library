#include "Matrix Visualization.h"
#include "matrix.h"
#include <iostream>
#include <complex>

namespace MatrixViz {
	//1.visualizeMultiplication() function implementation
	void visualizeMultiplication(const Matrix& mat1, const Matrix& mat2, const Matrix& result) {
		std::cout << "Multiplying matrices: \n\n";

		int maxRows = std::max({ mat1.getRows(), mat2.getCols(), result.getRows() });

		for (int i = 0; i < maxRows; ++i) {
			//Matrix A row
			if (i < mat1.getRows()) {
				std::cout << "|";
				for (int j = 0; j < mat1.getCols(); ++j) {
					std::cout << mat1.at(i, j) << " ";
				}
				std::cout << "|";
			}
			else {
				std::cout << std::string(mat1.getCols() * 2 + 3, ' ');
			}
			 
			// Multiplication sign
			if (i == mat1.getRows() / 2) {
				std::cout << " * ";
			}
			else {
				std::cout << "   ";
			}

			//Matrix B row
			if (i < mat2.getRows()) {
				std::cout << "|";
				for (int j = 0; j < mat2.getCols(); ++j) {
					std::cout << mat2.at(i, j) << " ";
				}
				std::cout << "|";
			}
			else {
				std::cout << std::string(mat2.getCols() * 2 + 3, ' ');
			}

			// Equals sign
			if (i == mat1.getRows() / 2) {
				std::cout << " = ";
			}
			else {
				std::cout << "   ";
			}

			//Result matrix row
			if (i < result.getRows()) {
				std::cout << "|";
				for (int j = 0; j < result.getCols(); ++j) {
					std::cout << result.at(i, j) << " ";
				}
				std::cout << "|";
			}
			else {
				std::cout << std::string(result.getCols() * 2 + 3, ' ');
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	//2.visualizeAddition() function implementation
	void visualizeAddition(const Matrix& mat1, const Matrix& mat2, const Matrix& result) {
		std::cout << "Adding matrices: \n\n";
		for (int i = 0; i < mat1.getRows(); ++i) {
			std::cout << "|";
			for (int j = 0; j < mat1.getCols(); ++j) {
				std::cout << mat1.at(i, j) << " ";
			}
			std::cout << "|";
			if (i == mat1.getRows() / 2) {
				std::cout << " + ";
			}
			else {
				std::cout << "   ";
			}
			std::cout << "|";
			for (int j = 0; j < mat2.getCols(); ++j) {
				std::cout << mat2.at(i, j) << " ";
			}
			std::cout << "|";
			if (i == mat1.getRows() / 2) {
				std::cout << " = ";
			}
			else {
				std::cout << "   ";
			}
			std::cout << "|";
			for (int j = 0; j < result.getCols(); ++j) {
				std::cout << result.at(i, j) << " ";
			}
			std::cout << "|" << std::endl;
		}
	}

	//3.visualizeSubtraction() function implementation
	void visualizeSubtraction(const Matrix& mat1, const Matrix& mat2, const Matrix& result) {
		std::cout << "Subtracting matrices:\n\n";
		for (int i = 0; i < mat1.getRows(); ++i) {
			std::cout << "|";
			for (int j = 0; j < mat1.getCols(); ++j) {
				std::cout << mat1.at(i, j) << " ";
			}
			std::cout << "|";
			if (i == mat1.getRows() / 2) {
				std::cout << " - ";
			}
			else {
				std::cout << "   ";
			}
			std::cout << "|";
			for (int j = 0; j < mat2.getCols(); ++j) {
				std::cout << mat2.at(i, j) << " ";
			}
			std::cout << "|";

			if (i == mat1.getRows() / 2) {
				std::cout << " = ";
			}
			else {
				std::cout << "   ";
			}

			std::cout << "|";
			for (int j = 0; j < result.getCols(); ++j) {
				std::cout << result.at(i, j) << " ";
			}
			std::cout << "|" << std::endl;

		}
		std::cout << std::endl;
	}
	

	//4.visualizeTranspose() function implementation
	void visualizeTranspose(const Matrix& mat, const Matrix& result) {
		std::cout << "Transposing matrix:" << std::endl;
		int maxRows = std::max(mat.getRows(), result.getRows());
		std::cout << "       A -----> A ^ T" << std::endl;
		std::cout << std::endl;

		for (int i = 0; i < maxRows; ++i) {

			//Print original matrix row
			if (i < mat.getRows()) {
				std::cout << "|";
				for (int j = 0; j < mat.getCols(); ++j) {
					std::cout << mat.at(i, j) << " ";
				}
				std::cout << "|";
			}
			else {
				//Blank space for allignment
				std::cout << "             ";
			}
			//Print arrow
			if (i == mat.getRows() / 2) {
				std::cout << " -----> ";
			}
			else {
				std::cout << "        ";
			}

			//Print transposed matrix row
			if (i < result.getRows()) {
				std::cout << "|";
				for (int j = 0; j < result.getCols(); ++j) {
					std::cout << result.at(i, j) << " ";
				}
				std::cout << "|";
			}
			std::cout << std::endl;
		}
	}

	//5.visualizeInverse() function implementation
	void visualizeInverse(const Matrix& mat, const Matrix& inverseMat){
		std::cout << "Calculating inverse of matrix:" << std::endl;
		std::cout << "Matrix:" << std::endl;
		mat.print();
		std::cout << "The inverse is: " << std::endl;
		inverseMat.print();
	}

	//6.visualizeEigenvalues() function implementation
	void visualizeEigenvalues(const Matrix& mat, const std::vector<std::complex<double>>& eigenvalues) {
		std::cout << "Calculating eigenvalues of matrix:" << std::endl;
		mat.print();

		if (mat.getRows() == 1) {
			std::cout << "Eigenvalue is the only element in the matrix" << std::endl;
		}
		if (mat.getRows() == 2) {
			std::cout << "Eigenvalues are the roots of the characteristic polynomial" << std::endl;
			std::cout << "Using quadratic formula y^2 - (" << mat.at(0,0)<<" + "<<mat.at(1,1)<<")y + ("<<mat.at(0,0)<<"*"<<mat.at(1,1)<<" - "<<mat.at(0,1)<<"*"<<mat.at(1,0)<<") = 0 to find eigenvalues where y represents the eigenvalues " << std::endl;
		}


		std::cout << "The eigenvalues are:" << std::endl;
		for (int i = 0; i < eigenvalues.size(); ++i) {
			if (eigenvalues[i].imag() != 0) {
				if (eigenvalues[i].imag() == 1)
					std::cout << eigenvalues[i].real() << " + i" << std::endl;
				else if (eigenvalues[i].imag() == -1)
					std::cout << eigenvalues[i].real() << " - i" << std::endl;
				else
				std::cout << eigenvalues[i].real() << " + " << eigenvalues[i].imag() << "i" << std::endl;
			}
			else {
				std::cout << eigenvalues[i].real() << std::endl;
			}
		}
	}

	//7.visualizeLUDecomposition() function implementation
	void visualizeLUDecomposition(const Matrix& mat, const Matrix& L, const Matrix& U) {
		std::cout << "LU Decomposition of matrix:" << std::endl;
		mat.print();
		std::cout << "L matrix:" << std::endl;
		L.print();
		std::cout << "U matrix:" << std::endl;
		U.print();
		std::cout << "Do you want to confirm L*U equals the original matrix? Press 1 for yes and 0 for no." << std::endl;
		int confirm;
		std::cin >> confirm;
		if (confirm == 1) {
			Matrix result = multiply(L, U);
			for (int i = 0; i < mat.getRows(); ++i) {
				for (int j = 0; j < mat.getCols(); ++j) {
					if (std::abs(result.at(i, j) - mat.at(i, j)) > 1e-10) {
						std::cout << "L*U does not equal the original matrix." << std::endl;

					}
					else {
						MatrixViz::visualizeMultiplication(L, U, result);
						std::cout << "L*U equals the original matrix." << std::endl;
					}
				}
			}
		}
		else if (confirm == 0) {
			std::cout << "You chose not to confirm L*U equals the original matrix." << std::endl;
		}
		else if (confirm != 0 && confirm != 1) {
			std::cout << "Invalid option. Please try again." << std::endl;
		}
	}

	//8.visualizeQRDecomposition() function implementation
	void visualizeQRDecomposition(const Matrix& mat, const Matrix& Q, const Matrix& R) {
		std::cout << "QR Decomposition of matrix:" << std::endl;
		mat.print();
		std::cout << "Q matrix:" << std::endl;
		Q.print();
		std::cout << "R matrix:" << std::endl;
		R.print();
		std::cout << "Do you want to confirm Q*R equals the original matrix? Press 1 for yes and 0 for no." << std::endl;
		int confirm;
		std::cin >> confirm;
		if (confirm == 1) {
			Matrix result = multiply(Q, R);
			for (int i = 0; i < mat.getRows(); ++i) {
				for (int j = 0; j < mat.getCols(); ++j) {
					if (std::abs(result.at(i, j) - mat.at(i, j)) > 1e-10) {
						std::cout << "Q*R does not equal the original matrix." << std::endl;
					}
					else {
						MatrixViz::visualizeMultiplication(Q, R, result);
					}
				}
			}
		}
		else if (confirm == 0) {
			std::cout << "You chose not to confirm Q*R equals the original matrix." << std::endl;
		}
		else if (confirm != 0 && confirm != 1) {
			std::cout << "Invalid option. Please try again." << std::endl;
		}
	}

	//9.visualizeCholeskyDecomposition() function implementation
	void visualizeCholeskyDecomposition(const Matrix& mat, const Matrix& L) {
		std::cout << "Cholesky Decomposition of matrix:" << std::endl;
		mat.print();
		std::cout << "L matrix:" << std::endl;
		L.print();
	}

	//10.visualizeLinearEquationSolution() function implementation
	void visualizeLinearEquationSolution(const Matrix& A, const std::vector<double>& x, const std::vector<double>& b) {
		std::pair<Matrix, Matrix> result = LUDecompose(A);
		Matrix L = result.first;
		Matrix U = result.second;
		std::cout << "Solving linear equations Ax = b. Remember we're solving for x" << std::endl;
		for (int i = 0; i < A.getRows(); ++i) {
			std::cout << "|";
			for (int j = 0; j < A.getCols(); ++j) {
				std::cout << "" << A.at(i, j) << " ";
			}
			std::cout << "|";
			if (i == A.getRows() / 2) {
				std::cout << " * ";
			}
			else {
				std::cout << "   ";
			}
			std::cout<< "|x|";
			if (i == A.getRows() / 2) {
				std::cout << " = ";
			}
			else {
				std::cout << "   ";
			}
			std::cout << "  |" << b[i] << "|" << std::endl;
		}
		std::cout << "LU decomposition is perfomed on matrix A and yields matrices L and U" << std::endl;
		std::cout << "Matrix L:" << std::endl;
		L.print();
		std::cout << "Matrix U:" << std::endl;
		U.print();
		std::cout << "Forward substitution is performed on L and b to get y. It solves the equation Ly = b" << std::endl;
		std::cout << "Vector b:" << std::endl;
		for (int i = 0; i < b.size(); ++i) {
			std::cout << b[i] << std::endl;
		}
		std::cout << "Backward substitution is performed on U and y to get x. It solves the equation Ux = y" << std::endl;
		std::cout << "The solution x in the order of the unknowns is:" << std::endl;
		for (int i = 0; i < x.size(); ++i) {
			std::cout << x[i] << std::endl;
		}
	}
}

//End of the source file