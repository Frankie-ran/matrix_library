#include<complex>
#include"matrix.h" //This is the header file for matrix operations
#include"Matrix Visualization.h" //This is the header file for matrix visualization functions
#include<iostream>
int main() {
	int option;
	do {
		std::cout << "Matrix operations menu:" << std::endl;
		std::cout << "Enter option. Press 0 to exit." << std::endl;
		std::cout << "1. Matrix addition" << std::endl;
		std::cout << "2. Matrix subtraction" << std::endl;
		std::cout << "3. Matrix multiplication" << std::endl;
		std::cout << "4. Matrix transpose" << std::endl;
		std::cout << "5. Matrix determinant" << std::endl;
		std::cout << "6. Matrix inverse" << std::endl;
		std::cout << "7. Eigenvalues" << std::endl;
		std::cout << "8. LU decomposition" << std::endl;
		std::cout << "9. QR decomposition" << std::endl;
		std::cout << "10. Cholesky decomposition" << std::endl;
		std::cout << "11. Solve linear equations" << std::endl;
		std::cout << "12.Check matrix types" << std::endl;
		std::cout << "13.Clear screen" << std::endl;

		std::cin >> option;
		switch (option) {
		case 0:
			std::cout << "Exiting program." << std::endl;
			break;

		case 1: {
			std::cout << "Matrix addition selected. The two matrices should have same dimensions" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows  for the matrices:" << std::endl;
			std::cin >> rows;
			std::cout << "Enter number of columns for the matrices:" << std::endl;
			std::cin >> cols;
			Matrix mat1(rows, cols);
			Matrix mat2(rows, cols);
			std::cout << "Enter elements of first matrix. You're entering the elements 1st row first" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat1.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat1.print();
			std::cout << "Enter elements of second matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat2.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat2.print();
			try {
				Matrix result = addMatrix(mat1, mat2);
				MatrixViz::visualizeAddition(mat1, mat2, result);
			}
			catch (const std::exception& e) {
				std::cerr << "Addition failed: " << e.what() << std::endl;
			}
			break;
		}
		case 2: {
			std::cout << "Matrix subtraction selected. The two matrices should have same dimensions" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows for the matrices:" << std::endl;
			std::cin >> rows;
			std::cout << "Enter number of columns for the matrices:" << std::endl;
			std::cin >> cols;
			Matrix mat1(rows, cols);
			Matrix mat2(rows, cols);
			std::cout << "Enter elements of first matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat1.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat1.print();
			std::cout << "Enter elements of second matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat2.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat2.print();
			try{
			Matrix result = subtraction(mat1, mat2);
			MatrixViz::visualizeSubtraction(mat1, mat2, result);
			}
			catch (const std::exception& e) {
				std::cerr << "Subtraction failed: " << e.what() << std::endl;
			}
			break;
		}
		case 3: {
			std::cout << "Matrix multiplication selected. The two matrices should have compatible dimensions" << std::endl;
			int rows1, cols1, rows2, cols2;
			std::cout << "Enter number of rows and columns for the first matrix:" << std::endl;
			std::cin >> rows1 >> cols1;
			Matrix mat1(rows1, cols1);
			std::cout << "Enter elements of first matrix:" << std::endl;
			for (int i = 0; i < rows1; ++i) {
				for (int j = 0; j < cols1; ++j) {
					std::cin >> mat1.at(i, j);
				}
			}
			std::cout << "You entered the following matrix. Please confirm." << std::endl;
			mat1.print();
			std::cout << "Enter number of rows and columns for the second matrix:" << std::endl;
			std::cin >> rows2 >> cols2;
			Matrix mat2(rows2, cols2);
			std::cout << "Enter elements of second matrix:" << std::endl;
			for (int i = 0; i < rows2; ++i) {
				for (int j = 0; j < cols2; ++j) {
					std::cin >> mat2.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat2.print();
			try{
			Matrix result = multiply(mat1, mat2);
			MatrixViz::visualizeMultiplication(mat1, mat2, result);
			}
			catch (const std::exception& e) {
				std::cerr << "Multiplication failed: " << e.what() << std::endl;
			}
			break;
		}
		case 4: {
			std::cout << "Matrix transpose selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			std::cout << "Performing transpose of the matrix:" << std::endl;
			Matrix result = transpose(mat);
			MatrixViz::visualizeTranspose(mat, result);
			break;
		}
		case 5: {
			std::cout << "Matrix determinant selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			try {
				std::cout << "Performing determinant of the matrix:" << std::endl;
				double result = determinant(mat, true);
			}
			catch (const std::exception& e) {
				std::cerr << "Determinant calculation failed: " << e.what() << std::endl;
			}

			break;
		} 
		case 6: {
			std::cout << "Matrix inverse selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			std::cout << "Performing inverse of the matrix:" << std::endl;
			try{
			Matrix result = inverse(mat, true);
			MatrixViz::visualizeInverse(mat, result);
			}
			catch (const std::exception& e) {
				std::cerr << "Inverse calculation failed: " << e.what() << std::endl;
			}
			break;
		}
		case 7: {
			std::cout << "Eigenvalues selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			std::cout << "Performing eigenvalue calculation of the matrix:" << std::endl;
			try{
		    std::vector<std::complex<double>> result = eigenvalues(mat);
			MatrixViz::visualizeEigenvalues(mat, result);
			}
			catch (const std::exception& e) {
				std::cerr << "Eigenvalue calculation failed: " << e.what() << std::endl;
			}
			break;
		}
		case 8: {
			std::cout << "LU decomposition selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			std::cout << "Performing LU decomposition of the matrix:" << std::endl;
			try {
				auto LU = LUDecompose(mat);
				Matrix L = LU.first;
				Matrix U = LU.second;
				MatrixViz::visualizeLUDecomposition(mat, LU.first, LU.second);
			}
			catch (const std::exception& e) {
				std::cerr << "LU decomposition failed: " << e.what() << std::endl;
			}
				break;
			
		}
		case 9: {
			std::cout << "QR decomposition selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			std::cout << "Performing QR decomposition of the matrix:" << std::endl;
			try {
				std::pair<Matrix, Matrix> result = QRDecompose(mat);
				MatrixViz::visualizeQRDecomposition(mat, result.first, result.second);
			}
			catch (const std::exception& e) {
				std::cerr << "QR decomposition failed: " << e.what() << std::endl;
			}
			break;
		}
		case 10: {
			std::cout << "Cholesky decomposition selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			std::cout << "Performing Cholesky decomposition of the matrix:" << std::endl;
			try{
			Matrix result = choleskyDecompose(mat);
			MatrixViz::visualizeCholeskyDecomposition(mat, result);
			}
			catch (const std::exception& e) {
				std::cerr << "Cholesky decomposition failed: " << e.what() << std::endl;
			}
			break;
		}
		case 11: {
			std::cout << "Solving linear equations selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			std::vector<double> b(rows);
			std::cout << "Enter elements of vector b:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				std::cin >> b[i];
			}
			std::cout << "You entered the following vector:" << std::endl;
			for (int i = 0; i < b.size(); i++) {
				std::cout << b[i] << std::endl;
			}
			std::cout << std::endl;
			try{
			std::vector<double> result = solveLinearEquations(mat, b);
			MatrixViz::visualizeLinearEquationSolution(mat, result, b);
			}
			catch (const std::exception& e) {
				std::cerr << "Linear equation solving failed: " << e.what() << std::endl;
			}
			break;
		}
		case 12: {
			std::cout << "Check matrix types selected" << std::endl;
			int rows, cols;
			std::cout << "Enter number of rows and columns for the matrix:" << std::endl;
			std::cin >> rows >> cols;
			Matrix mat(rows, cols);
			std::cout << "Enter elements of the matrix:" << std::endl;
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					std::cin >> mat.at(i, j);
				}
			}
			std::cout << "You entered the following matrix:" << std::endl;
			mat.print();
			int checkingOption;
			do {
				std::cout << "Enter option. Press 0 to exit." << std::endl;
				std::cout << "1. Check if the matrix is symmetric" << std::endl;
				std::cout << "2. Check if the matrix is an identity matrix" << std::endl;
				std::cout << "3. Check if the matrix is a square matrix" << std::endl;
				std::cout << "4. Check if the matrix is a skew-symmetric matrix" << std::endl;
				std::cout << "5. Check if the matrix is a diagonal matrix" << std::endl;
				std::cout << "6. Check if the matrix is a upper triangular matrix" << std::endl;
				std::cout << "7. Check if the matrix is a lower triangular matrix" << std::endl;
				std::cout << "8. Check if the matrix is a orthogonal matrix" << std::endl;
				std::cout << "9. Check if the matrix is a singular matrix" << std::endl;
				std::cin >> checkingOption;
				switch (checkingOption) {
				case 0:
					std::cout << "Exiting program." << std::endl;
					break;
				case 1: {
					std::cout << "Checking if the matrix is symmetric:" << std::endl;
					if (isSymmetric(mat)) {
						std::cout << "The matrix is symmetric." << std::endl;
					}
					else {
						std::cout << "The matrix is not symmetric." << std::endl;
					}
					break;
				}
				case 2: {
					std::cout << "Checking if the matrix is an identity matrix:" << std::endl;
					if (isIdentity(mat)) {
						std::cout << "The matrix is an identity matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not an identity matrix." << std::endl;
					}
					break;
				}
				case 3: {
					std::cout << "Checking if the matrix is a square matrix:" << std::endl;
					if (isSquare(mat)) {
						std::cout << "The matrix is a square matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not a square matrix." << std::endl;
					}
					break;
				}
				case 4: {
					std::cout << "Checking if the matrix is a skew-symmetric matrix:" << std::endl;
					if (isSkewSymmetric(mat)) {
						std::cout << "The matrix is a skew-symmetric matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not a skew-symmetric matrix." << std::endl;
					}
					break;
				}
				case 5: {
					std::cout << "Checking if the matrix is a diagonal matrix:" << std::endl;
					if (isDiagonal(mat)) {
						std::cout << "The matrix is a diagonal matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not a diagonal matrix." << std::endl;
					}
					break;
				}
				case 6: {
					std::cout << "Checking if the matrix is a upper triangular matrix:" << std::endl;
					if (isUpperTriangular(mat)) {
						std::cout << "The matrix is an upper triangular matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not an upper triangular matrix." << std::endl;
					}
					break;
				}
				case 7: {
					std::cout << "Checking if the matrix is a lower triangular matrix:" << std::endl;
					if (isLowerTriangular(mat)) {
						std::cout << "The matrix is a lower triangular matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not a lower triangular matrix." << std::endl;
					}
					break;
				}
				case 8: {
					std::cout << "Checking if the matrix is a orthogonal matrix:" << std::endl;
					if (isOrthogonal(mat)) {
						std::cout << "The matrix is an orthogonal matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not an orthogonal matrix." << std::endl;
					}
					break;
				}
				case 9: {
					std::cout << "Checking if the matrix is a singular matrix:" << std::endl;
					if (isSingular(mat)) {
						std::cout << "The matrix is a singular matrix." << std::endl;
					}
					else {
						std::cout << "The matrix is not a singular matrix." << std::endl;
					}
				}
					  break;
				default:
					std::cout << "Invalid option. Please try again." << std::endl;
					break;
				}
			} while (checkingOption != 0);
			break;
		}
		case 13: {
			std::cout << "Clearing screen." << std::endl;
			system("cls");
			break;
		}

		case 14:{
			Matrix mat(2, 2);

			mat.print();
		}
		default:
			std::cout << "Invalid option. Please try again." << std::endl;
			break;
		}
	} while (option != 0);
	return 0;
}