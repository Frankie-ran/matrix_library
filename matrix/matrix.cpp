#include"matrix.h"
#include<iostream>
#include<complex>
#include<cmath>
#include<utility>
#include<vector>
#include<stdexcept>
#include "Matrix Visualization.h"


//Matrix constructor implementation and initialization of all values to 0.0
Matrix::Matrix(int r, int c) {
	rows = r;
	cols = c;
	data = new double* [rows];
	for (int i = 0; i < rows; i++) {
		data[i] = new double[cols];
		for (int j = 0; j < cols; j++) {
			data[i][j] = 0.0;
		}
	}
}

//Identity matrix constructor implementation
Matrix identityMatrix(int size) {
	Matrix identity(size, size); //Create a square matrix of size x size
	for (int i = 0; i < size; ++i) {
		identity.at(i, i) = 1.0; //Set diagonal elements to 1
	}
	return identity; //Return the identity matrix
}

//Copy constructor implementation
Matrix::Matrix(const Matrix& other) {
	rows = other.rows;
	cols = other.cols;
	data = new double* [rows];
	for (int i = 0; i < rows; ++i) {
		data[i] = new double[cols];
		for (int j = 0; j < cols; ++j) {
			data[i][j] = other.data[i][j]; //Copy the values from the other matrix
		}
	}
}

//Read-only access to the matrix element at row r and column c
const double& Matrix::at(int r, int c) const {
	return data[r][c];
}


//Read and write access to the matrix element at row r and column c
double& Matrix::at(int r, int c) {
	return data[r][c];
}


//Returns the number of rows
int Matrix::getRows() const {
	return rows;
}


//Returns the number of columns
int Matrix::getCols() const {
	return cols;
}

//Print function implementation
void Matrix::print() const{
	for (int i = 0; i < getRows(); ++i) {
		std::cout << "|";
		for (int j = 0; j < getCols(); ++j) {
			std::cout <<""<< at(i, j) << " ";
		}
		std::cout << "|" << std::endl;
	}
}

//Destructor implementation
Matrix::~Matrix() {
	for (int i = 0; i < rows; i++) {
		delete[] data[i];
	}
	delete[] data;
}


//Matrix addition implementation
Matrix addMatrix(const Matrix& mat1, const Matrix& mat2) {
	//Check if dimensions match
	if (mat1.getRows() != mat2.getRows() || mat1.getCols() != mat2.getCols()) {
		throw
			std::invalid_argument("Matrices dimensions must match for addition");
	}
	else {
		//Create a result matrix of the same size
		Matrix resultMatrix(mat1.getRows(), mat1.getCols());

		// Perform the addition of corresponding elements
		for (int i = 0; i < mat1.getRows(); i++) {
			for (int j = 0; j < mat1.getCols(); j++) {
				resultMatrix.at(i, j) = mat1.at(i, j) + mat2.at(i, j);
			}
		}
		return resultMatrix;
	}
}


//Matrix subtraction implementation
Matrix subtraction(const Matrix& mat1, const Matrix& mat2) {
	//Check if dimensions match
	if (mat1.getRows() != mat2.getRows() || mat1.getCols() != mat2.getCols()) {
		throw
			std::invalid_argument("Matrices dimensions must match for addition");
	}
	else {
		//Create a result matrix of the same size
		Matrix resultMatrix(mat1.getRows(), mat1.getCols());

		// Perform the subtraction of corresponding elements
		for (int i = 0; i < mat1.getRows(); i++) {
			for (int j = 0; j < mat1.getCols(); j++) {
				resultMatrix.at(i, j) = mat1.at(i, j) - mat2.at(i, j);
			}
		}
		return resultMatrix;
	}
}


//Matrix multiplication implementation
Matrix multiply(const Matrix& mat1, const Matrix& mat2) {
	//Check whether number ofcolumns  in mat1 is equal to number of rows in mat2
	if(mat1.getCols() != mat2.getRows()){
		throw
			std::invalid_argument("Matrix multiplication not possible. Incompatible dimensions");
	}else{
		//Create a result matrix
		Matrix resultMatrix(mat1.getRows(), mat2.getCols());

		//Perform matrix multiplication
		for (int i = 0; i < mat1.getRows(); ++i) {
			for (int j = 0; j < mat2.getCols(); ++j) {
				double sum = 0.0;
				for (int k = 0; k < mat1.getCols(); ++k) {
					sum += mat1.at(i, k) * mat2.at(k, j);
				} 
				resultMatrix.at(i, j) = sum;
			}
		}
		return resultMatrix;
	}

}


//Matrix transpose implementation
Matrix transpose(const Matrix& mat)  { 
	Matrix resultMatrix(mat.getCols(), mat.getRows());
	for (int i = 0; i < mat.getRows(); i++) {
		for (int j = 0; j < mat.getCols(); j++) {
			resultMatrix.at(i, j) = mat.at(j, i);
		}
}
	return resultMatrix;
}


//GetMinor function to use in the determinant and cofactor functions
static Matrix getMinor(const Matrix& mat, int rowToRemove, int colToRemove ) {
	int rows = mat.getRows();
	int cols = mat.getCols();

	//Construct the minor matrix
	Matrix minor(rows - 1, cols - 1); //Create a minor matrix of size (rows-1) x (cols-1)
	int minorRow = 0;

	for (int i = 0; i < rows; ++i) {
		if (i == rowToRemove)
			continue; //skip the row
		int minorCol = 0;
		for (int j = 0; j < cols; ++j) {
			if (j == colToRemove)
				continue; //Skip the column
			minor.at(minorRow, minorCol) = mat.at(i, j);
			++minorCol;
		}
		++minorRow;
	}
	return minor;
}


//Matrix determinant implementation
//Also used in cofactor function implementation
double determinant(const Matrix& mat, bool visualize, int depth) {
	int n = mat.getRows();
	if (n != mat.getCols()) {
		throw
			std::invalid_argument("NUmber of rows must be equal to number of columns");
	}

	// Base case for 1x1  matrices
		if (n == 1) {
			if (visualize) {
				std::cout <<std::string(depth*2, ' ')<< "Determinant of 1x1 matrix: " << mat.at(0, 0) << std::endl;
			}
			return mat.at(0, 0);
		}

		// Base case for 2x2 matrices
	    if (n == 2) {
			double det = mat.at(0, 0) * mat.at(1, 1) - mat.at(0, 1) * mat.at(1, 0);
			if (visualize) {
				std::cout << std::string(depth * 2, ' ')
					<< "Base case for 2x2 matrix: (" << mat.at(0, 0) 
					<< " * " << mat.at(1, 1) << ") - (" << mat.at(0, 1) 
					<< " * " << mat.at(1, 0) << ") = " << det << std::endl;
			}
			return det;
		}

		// Recursive case for larger matrices
			double det = 0.0;
			for (int j = 0; j < n; ++j) {
				Matrix minor = getMinor(mat, 0, j);

				if (visualize) {
					std::cout<< std::string(depth*2, ' ') <<
						"Expanding along row 0, column " << j << mat.at(0, j) << "\n";
					std::cout << std::string(depth * 2, ' ') << "Minor matrix at depth:" << depth << "\n";
					minor.print();
				}
				double minorDet = determinant(minor, visualize, depth + 1);

				double cofactor = ((j % 2 == 0) ? 1 : -1) * mat.at(0, j) *
					determinant(minor, false);
				det += cofactor;

				if (visualize) {
					std::cout<<std::string(depth * 2, ' ')
						<< "Cofactor contribution: " << cofactor << "\n";

				}
			}
			if (visualize && depth == 0) {
				std::cout << std::string(depth * 2, ' ')
					<< "Determinant of the matrix: " << det << std::endl;
			}
			return det;	
}

//Matrix cofactor function to use in the inverse function
static Matrix cofactorMatrix(const Matrix& mat){
	//Square matrix validated at the inverse function
	int n = mat.getRows();

	Matrix cofactors(n, n); //Create a cofactor matrix.

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			Matrix minor = getMinor(mat, i, j);
			double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
			cofactors.at(i, j) = sign * determinant(minor, false);
		}
	}
	return cofactors;

}

//Matrix inverse implementation
Matrix inverse(const Matrix& mat, bool visualize) {
	int n = mat.getRows();
	if (n != mat.getCols()) {
		throw
			std::invalid_argument("Only square matrices can be inverted");
	}
	else {
		double det = determinant(mat, false);
		if (det == 0.0) {
			throw
				std::runtime_error("Matrix is not invertible. Its determinant is zero");
		}
		else {
			//Create a cofactor matrix
			Matrix cofactors = cofactorMatrix(mat);
			if (visualize) {
				std::cout << "Cofactor matrix:" << std::endl;
				cofactors.print();
			}


			//Create an adjugate of the cofactor matrix
			Matrix adjugate = transpose(cofactors);
			if (visualize) {
				MatrixViz::visualizeTranspose(cofactors, adjugate);
			}

			//Create the inverse matrix
			Matrix inverseMat(n, n);

			//Divide each element of the matrix by the determinant
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					inverseMat.at(i, j) = adjugate.at(i, j) / det;
				}
			}
			if (visualize) {
				std::cout << "Inverse matrix:" << std::endl;
				inverseMat.print();
			}
			return inverseMat;
		}
	}
}

//Eigenvalues function implementation
std::vector<std::complex<double>>
	eigenvalues(const Matrix& mat) {
	if (mat.getRows() != mat.getCols()) {
		throw
			std::invalid_argument("Eigenvalues calculation can only be done on square matrices");
	}
	if(mat.getRows() == 1 && mat.getCols() == 1) { 
		//This handles matrices of 1 by 1 size.
		//Eigenvalue is the only element in the matrix
		double value = mat.at(0, 0);
		return
		{ std::complex<double>(value,0),
		std::complex<double>(value,0) };
	
	}
	if (mat.getRows() == 2 && mat.getCols() == 2) {
		//This handles matrices of 2 by 2 size
		//Uses quadratic formula to find eigenvalues
		double a = mat.at(0, 0);
		double b = mat.at(0, 1);
		double c = mat.at(1, 0);
		double d = mat.at(1, 1);

		//Calculate discriminant and trace
		double trace = a + d;
		double determinant = a * d - b * c;
		double discriminant = trace * trace - 4 * determinant;

		std:: complex<double>
			sqrt_discriminant =
			std::sqrt(std::complex<double>(discriminant, 0));

		std::complex<double> lambda1 = (trace + sqrt_discriminant) / 2.0;
		std::complex<double> lambda2 = (trace - sqrt_discriminant) / 2.0;

		return{ lambda1, lambda2 };
	}
	 

		//This handles matrices of size greater than 2 by 2
		//It uses QRDecompose() repeatedly to converge to a matrix whose diagonal elements are the eigenvalues
		const int maxIterations = 1000;
		const double tolerance = 1e-10;

		Matrix A = mat;
		for (int iter = 0; iter < maxIterations; ++iter) {
			std::pair<Matrix, Matrix> result = QRDecompose(A);
			Matrix Q = result.first;
			Matrix R = result.second;
			A = multiply(R, Q); //A = RQ

			//Check for convergence
			bool converged = true;
			for (int i = 0; i < A.getRows(); ++i) {
				for (int j = 0; j < A.getCols(); ++j) {
					if (i != j && std::abs(A.at(i, j)) > tolerance) {
						converged = false;
						break;
					}
				}
				if (!converged) break;
				}
			if (converged) break;
         }
		//Extract eigenvalues from the diagonal of the converged matrix
		std::vector<std::complex<double>> result;
		for (int i = 0; i < A.getRows(); ++i) {
			result.push_back(A.at(i, i));
		}
		return result;

	
}

//LU decomposition
std::pair<Matrix, Matrix>
LUDecompose(const Matrix& mat) {
	int n = mat.getRows();
	if (n != mat.getCols()) {
		throw
			std::invalid_argument("LU decomposition requires a square matrix");
	}
	for (int i = 0; i < n; ++i) {
		if (std::abs(mat.at(i, i)) < 1e-10) {
			throw
				std::runtime_error("Zero pivot encountered. Pivoting hasn't been implemented so do not input a matrix with zero at the diagonal");
		}
	}

	Matrix L(n,n); //Create the lower matrix
	Matrix U(n,n); // Create the upper matrix

	for (int i = 0; i < n; ++i) {
		//Compute row i of U using the formula U(i,j) = A(i,j) - sum(L(i,k)*U(k,j)) where k is from 0 to i-1
		for (int j = i; j < n; ++j) {
			double sum = 0.0;
			for (int k = 0; k < i; ++k) {
				sum += L.at(i, k) * U.at(k, j);
			}
			U.at(i, j) = mat.at(i, j) - sum;
		}

		L.at(i, i) = 1.0; // Set the diagonal of L to 1.0

		//Compute column i of L using the formula L(j,i) = (A(j,i) - sum(L(j,k)*U(k,i))) / U(i,i) where k is from 0 to i-1
		for (int j = i + 1; j < n; ++j) {
			double sum = 0.0;
			for (int k = 0; k < i; ++k) {
				sum += L.at(j, k) * U.at(k, i);
			}
			if (U.at(i, i) == 0.0) {
				throw
					std::runtime_error("Zero pivot encountered at U");
			}
			L.at(j, i) = (mat.at(j, i) - sum) / U.at(i, i);
		}
	}
	return std::make_pair(L, U);
	//return{L,U}; This is another way to return two values
}

//2. QR decomposition
std::pair<Matrix, Matrix>
QRDecompose(const Matrix& mat) {
	int rows = mat.getRows();
	int cols = mat.getCols();
	if (rows < cols) {
		throw
			std::invalid_argument("QR decomposition requires the number of rows to be equal or more than the number of columns");
	}


	 //Create Orthonormal Matrix Q
	Matrix Q(rows, cols);
	//Create upper triangular matrix R
	Matrix R(cols, cols);

	for (int j = 0; j < cols; ++j) {
		//Copy column j of mat to q
		for (int i = 0; i < rows; ++i) {
			Q.at(i, j) = mat.at(i, j);
		}
		//Orhorgogonalize against previous columns
		for (int k = 0; k < j; ++k) {
			double dot_product = 0.0;
			for (int i = 0; i < rows; ++i) {
				dot_product += mat.at(i, j) * Q.at(i, k);//This will be R.at(k, j)
			}
			R.at(k, j) = dot_product;
			for (int i = 0; i < rows; ++i) {
				Q.at(i, j) -= dot_product * Q.at(i, k);
			}
		}
		//Normalize the column
		double norm = 0.0;
		for (int i = 0; i < rows; ++i) {
			norm += Q.at(i, j) * Q.at(i, j);
		}
		norm = std::sqrt(norm);
		if (norm == 0.0) {
			throw
				std::runtime_error("Matrix contains linearly dependent columns");
		}
		R.at(j, j) = norm;
		for (int i = 0; i < rows; ++i) {
			Q.at(i, j) /= norm;
		}

	}
	return std::make_pair(Q, R);
}


//A helper function to check whether a matrix is eligeble for cholesky decomposition
static bool isCholeskyEligible(const Matrix& mat) {
	int n = mat.getRows();
	if (n != mat.getCols()) {
		return false; //Not a square matrix
	}
	//Check if the matrix is symmetric
	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			if (mat.at(i, j) != mat.at(j, i)) {
				return false; //Not symmetric
			}
		}
	}
	//Check if eigenvalues are real and positive
	auto eigen = eigenvalues(mat);
	for (const auto& val : eigen) {
		if (val.real() <= 0) {
			return false; //Not positive definite
		}
		if (std::abs(val.imag()) > 1e-10) {
			return false; //Not real
		}
	}
	return true; //Matrix is eligible for Cholesky decomposition
}



//Cholesky decomposition
 Matrix choleskyDecompose(const Matrix& mat) {
	if (!isCholeskyEligible(mat)) {
		throw
			std::invalid_argument("Matrix is not eligible for Cholesky decomposition");
	}
	int n = mat.getRows();
	//Create the result matrix L
	Matrix L(n, n); //Lower triangular matrix


	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <= i; ++j) {
			double sum = 0.0;
			for (int k = 0; k < j; ++k) {
				sum += L.at(i, k) * L.at(j, k);
			}
			if (i == j) {
				L.at(i, j) = std::sqrt(mat.at(i, i) - sum);
			}
			else {
				L.at(i, j) = (mat.at(i, j) - sum) / L.at(j, j);
			}
		}
	}
	return L;
}


 //Linear equation solving functionality
 std::vector<double>
	 solveLinearEquations(const Matrix& A, const std::vector<double>& b) {
	 int r = A.getRows();
	 int c = A.getCols();
	 //Check if A is square and the number of rows in A matches the size of b
	 if (r != c) {
		 throw
			 std::invalid_argument("Matrix A must be square");
	 }
	 if (b.size() != r) {
		 throw
			 std::invalid_argument("Matrix A and vector b must have compatible dimensions");
	 }
	 //Step 1: LU decomposition
	 std::pair<Matrix, Matrix>
		 LU = LUDecompose(A);
	 Matrix L = LU.first;
	 Matrix U = LU.second;
     //Create the solution vectors
	 std::vector<double> y(r, 0.0);//Creates a vector of size r to store intermediate results of forward substitution when solving Ly = b
	 std::vector<double> x(r, 0.0);//Creates a vector of size r to store the final solution when solving Ux = y

	 //Step 2 : Forward substitution to solve Ly = b
	 for (int i = 0; i < r; ++i) {
		 y[i] = b[i];
		 for (int j = 0; j < i; ++j) {
			 y[i] -= L.at(i, j) * y[j];
		 }
		 y[i] /= L.at(i, i);
	 }

	 //Step 3: Backward substitution to solve Ux = y
	 for (int i = r - 1; i >= 0; --i) {
		 x[i] = y[i];
		 for (int j = i + 1; j < r; ++j) {
			 x[i] -= U.at(i, j) * x[j];
		 }
		 x[i] /= U.at(i, i);
	 }
	 return x; //Return the solution 
 }



 //Different matrix implementations:
 //1. Identity matrix
 bool isIdentity(const Matrix& mat) {
	 double tolerance = 1e-9;
	 int r = mat.getRows();
	 int c = mat.getCols();
	 //Check if the matrix is square
	 if (r != c) {
		 return false; //Not a square matrix so not an identity matrix
	 }
	 else {
		 for (int i = 0; i < r; ++i) {
			 for (int j = 0; j < c; ++j) {
				 if (i == j && std::abs(mat.at(i, j) - 1.0) > tolerance) {
					 return false;
				 }
				 else if(i != j && std::abs(mat.at(i, j)) > tolerance) {
					 return false;
				 }
			 }
		 }
	 }
	 return true;
 }

 //2. Symmetric matrix
 bool isSymmetric(const Matrix& mat) {
	 if (mat.getRows() != mat.getCols()) {
		 return false; //Not a square matrix so cannot be symmetric
	 }
	 else {
		 for (int i = 0; i < mat.getRows(); ++i) {
			 for (int j = 0; j < mat.getCols(); ++j) {
				 if (mat.at(i, j) != mat.at(j, i)) {
					 return false; //Not symmetric
				 }
			 }
		 }
	 }
	 return true;
 }

 //3.Square matrix
 bool isSquare(const Matrix& mat) {
	 return mat.getRows() == mat.getCols();
 }

 
 //4. Skew-symmetric matrix
 bool isSkewSymmetric(const Matrix& mat) {
	 if (mat.getRows() != mat.getCols()) {
		 return false; //Not a square matrix so cannot be skew-symmetric
	 }
	 else {
		 for (int i = 0; i < mat.getRows(); ++i) {
			 for (int j = 0; j < mat.getCols(); ++j) {
				 if (mat.at(i, j) != mat.at(j, i) / -1) {
					 return false; //Not skew-symmetric
				 }
			 }
		 }
	 }
	 return true;
 }

 //5.Diagonal matrix
 bool isDiagonal(const Matrix& mat) {
	 if (mat.getRows() != mat.getCols()) {
		 return false; //Not square. Cannot be diagonal
	 }
	 else {
		 for (int i = 0; i < mat.getRows(); ++i) {
			 for (int j = 0; j < mat.getCols(); ++j) {
				 if (i != j && mat.at(i, j) != 0.0) {
					 return false; //Off-diagonal elements must be zero
				 }
			 }
		 }
	 }
	 return true;
 }

 //6. Upper triangular matrix
 bool isUpperTriangular(const Matrix& mat) {
	 int r = mat.getRows();
	 int c = mat.getCols();
	 if (r != c) {
		 return false; //Not square. Cannot be upper triangular
	 }
	 else {
		 for (int i = 1; i < r; ++i) {
			 for (int j = 0; j < i; ++j) {
				 if (mat.at(i, j) != 0.0) {
					 return false; //Elements below the diagonal must be zero
				 }
			 }
		 }
	 }
	 return true; 
 }

 //7. Lower triangular matrix
 bool isLowerTriangular(const Matrix& mat) {
	 int r = mat.getRows();
	 int c = mat.getCols();
	 if (r != c) {
		 return false; //Not square. Cannot be lower triangular
	 }
	 else {
		 for (int i = 0; i < r; ++i) {
			 for (int j = i + 1; j < c; ++j) {
				 if (mat.at(i, j) != 0.0) {
					 return false; //Elements above the diagonal must be zero
				 }
			 }
		 }
	 }
	 return true;
 }

 //8. Orthogonal matrix
 bool isOrthogonal(const Matrix& mat) {
	 int r = mat.getRows();
	 int c = mat.getCols();
	 if (r != c) {
		 return false; //Not square. Cannot be orthogonal
	 }
	 else {
		 Matrix transposeMat = transpose(mat);
		 Matrix identity = multiply(mat, transposeMat);
		 return isIdentity(identity); //Check if the product is an identity matrix
	 }
 }

 bool isSingular(const Matrix& mat) {
	 //Here we're checking to see if the diagonal of U after LU decomposition has zero. If it does, the matrix is singular
	 if (mat.getRows() != mat.getCols()) {
		 throw
			 std::invalid_argument("Matrix singularity is only defined in square matrices");
	 }
	 double tolerance = 1e-9;
	 std::pair<Matrix, Matrix>
	  LU = LUDecompose(mat);
	 Matrix L = LU.first;
	 Matrix U = LU.second;

	 int r = mat.getRows();
	 for (int i = 0; i < r; ++i) {
		 if (std::abs(U.at(i, i)) < tolerance) {
			 return true;
		 }
	 }
	 return false;
 }


 //End of  source file