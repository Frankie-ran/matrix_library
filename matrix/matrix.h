#pragma once
#include <iostream>
#include <complex> //This library is used since eigenvalues might be complex numbers.
#include<utility> // For std::pair. This library is used when a function rerurns two values.
#include <cmath> //For sqrt function
#include <vector> //For std::vector. This library is used to return a vector of eigenvalues
class Matrix {
private:
	int rows, cols; // Number of rows and columns
	double** data; // POinter to a 2D array
public:
	//Matrix cconstructor declaration
	Matrix(int r, int c);// Constructor to create a matrix of size r x c. Initializes all values to 0.0

	
	//Access functions
	double& at(int r, int c); // Returns a reference to the matrix element at row r and column c. Allows modification of the matrix
	const double& at(int r, int c) const; //Used when matrix object is constant. Read only access

	//Matrix size functions
	int getRows() const; // Returns the number of rows
	int getCols() const; //Returns the number of columns

	//Matrix print function
	void print() const;

	Matrix(const Matrix& other); // Copy constructor. Initializes a new matrix as a copy of another matrix	

	//Destructor declaration
	~Matrix();
};
//Identity matrix constructor declaration
Matrix identityMatrix(int size); //Initializes a square matrix of size x size with 1s on the diagonal and zeros elsewhere


	//Simple matrix operations
	Matrix addMatrix(const Matrix& mat1, const Matrix& mat2); //Adding two matrices together. Both provided by the user
	Matrix multiply(const Matrix& mat1, const Matrix& mat2); //Multiply two matrices. Both provided by the user
	//Matrix multiplication is not commutative. So the order of the matrices matters. You're multiplying mat1 by mat2
	Matrix subtraction(const Matrix& mat1, const Matrix& mat2); // Subtract a matrix from another.You're subtracting mat2 from mat1. Both provided by the user
	Matrix transpose(const Matrix& mat); //Transpose a matrix provided by the user

	//Other Matrix operations
	double determinant(const Matrix& mat, bool visualize, int depth = 0); //Calculating a determinant of a matrix provided by the library user
	//If visualize is true, the function will print the steps of the calculation


	Matrix inverse(const Matrix& mat, bool visualize); //Finding the inverse of a matrix provided by the user

	//Function declaration for eigenvalues
	std::vector<std::complex<double>>
		eigenvalues(const Matrix& mat);

	//MATRIX DECOMPOSITION ALGORITHMS
	
	
	//1. LU decomposition. This decomposition returns matrix L and U
	//Matrix L(lower triangular matrix) has ones in its diagonal and possibly non-zero values below the diagonal. It has zeros above the diagonal.
	//Matrix U(Upper triangular matrix) has non zero values on the diagonal and above, with zeros below the diagonal
	std::pair<Matrix, Matrix>
		LUDecompose(const Matrix& mat);

	//2. QR decomposition. This decomposition returns matrices Q and R
	// Matrix Q. This matrix has its rows(or columns, depending on context) perpendicular to each other.
	// In the case of real numbers, the columns are also unit vectors
	//Matrix R. This matrix has all entries below the diagonal equal to zero
	std::pair<Matrix, Matrix>
		QRDecompose(const Matrix& mat);

	//3. Cholesky decomposition. This decomposition returns matrix L. 
	//Matrix L. This matrix is lower triangular and has positive values on the diagonal
	Matrix choleskyDecompose(const Matrix& mat);


//Linear equation solving functionality
	std::vector<double>
		solveLinearEquations(const Matrix& A, const std::vector<double>& b); //Solving a system of linear equations Ax = b
	//This function uses LU decomposition to solve the system of equations
	//It returns the solution vector x



//Matrix types implementation.
	bool isSymmetric(const Matrix& mat);//A symmetric matrix is the same as its transpose.
	bool isSquare(const Matrix& mat);//This matrix has the same number of rows as columns.
	bool isSkewSymmetric(const Matrix & mat);//A skew-symmetric matrix is the same as its negative transpose.
	bool isDiagonal(const Matrix& mat);//A diagonal matrix has non-zero values only on the diagonal.
	bool isUpperTriangular(const Matrix& mat);//An upper triangular is a square matrix that has non-zero values only on the diagonal and above.
	bool isLowerTriangular(const Matrix& mat);//A lower triangular matrix is a square matrix that has non-zero values only on the diagonal and below.
	bool isOrthogonal(const Matrix& mat);//An orthogonal matrix has its rows(or columns) perpendicular to each other. 
	//Its transpose is equal to its inverse.

	bool isIdentity(const Matrix& mat);//An identity matrix is a square matrix that has ones on the diagonal and zeros elsewhere. 
	//It has to be square too. 

	bool isSingular(const Matrix& mat);//A singular matrix is a square matrix that has a determinant of zero. It's not invertible.

	//End of the header file
