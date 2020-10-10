#include "Matrix.h"
#include "Functions.h"
#include <utility>
#include <iomanip>
#define PRECISION 6

/**
 * Constructor with allocating memory for matrix
 */
Matrix::Matrix(int n) {
	size = n;
	arr = new double* [size];
	for (int i = 0; i < size; ++i) {
		arr[i] = new double[size]{0};
	}
}

/**
 * Copy constructor
 */
Matrix::Matrix(const Matrix& other) : Matrix(other.size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			arr[i][j] = other.arr[i][j];
		}
	}
}

/**
 * Copy assignment operator
 */
Matrix& Matrix::operator=(Matrix other) {
	swap(*this, other);

	return *this;
}

/**
 * Move constructor
 */
Matrix::Matrix(Matrix&& other) noexcept {
	swap(*this, other);
}

/**
 * Destructor
 */
Matrix::~Matrix() {
	for (int i = 0; i < size; ++i) {
		delete[] arr[i];
	}

	delete[] arr;
}

/**
 * Swap
 */
void swap(Matrix& lhs, Matrix& rhs) noexcept {
	using std::swap;

	swap(lhs.size, rhs.size);
	swap(lhs.arr, rhs.arr);
}

/**
 * Input matrix from "fin" stream
 */
void Matrix::input(std::istream& fin) const {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			fin >> arr[i][j];
		}
	}
}

/**
 * Print matrix into "ostr" stream
 */
void Matrix::print(std::ostream& ostr) const {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			ostr << arr[i][j] << " ";
		}

		ostr << std::endl;
	}
}

/**
 * Indexer
 */
double*& Matrix::operator[](const int index) const {
	return arr[index];
}

/**
 * Operator for multiply matrix: like A*B
 */
Matrix Matrix::operator*(const Matrix& other) const {
	Matrix new_matrix(size);
	for (int i = 0; i < size; ++i) {
		for (int k = 0; k < size; ++k) {
			for (int j = 0; j < size; ++j) {
				new_matrix[i][k] += (*this)[i][j] * other[j][k];
			}
		}
	}

	return new_matrix;
}

/**
 * Operator for multiply matrix to vector: like A*x
 */
Vector Matrix::operator*(const Vector& vec) const {
	Vector result(size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			result[i] += (*this)[i][j] * vec[j];
		}
	}

	return result;
}

/**
 * First norm of matrix
 */
double Matrix::get_first_norm() const {
	double result = -1;
	for (int i = 0; i < size; ++i) {
		double sum = 0;
		for (int j = 0; j < size; ++j) {
			sum += fabs((*this)[i][j]);
		}

		if (sum > result) {
			result = sum;
		}
	}

	return result;
}

/**
 * Second norm of matrix
 */
double Matrix::get_second_norm() const {
	double result = -1;
	for (int i = 0; i < size; ++i) {
		double sum = 0;
		for (int j = 0; j < size; ++j) {
			sum += fabs((*this)[j][i]);
		}

		if (sum > result) {
			result = sum;
		}
	}

	return result;
}

/**
 * Third norm of matrix
 */
double Matrix::get_third_norm() const {
	Matrix tmp_matrix(size);
	// A * A^T
	for (int i = 0; i < size; ++i) {
		for (int k = 0; k < size; ++k) {
			for (int j = 0; j < size; ++j) {
				tmp_matrix[i][k] += (*this)[j][i] * (*this)[j][k];
			}
		}
	}

	// Diagonal matrix for A * A^T
	Matrix diagonal_matrix = tmp_matrix.get_diagonal_matrix();
	// Finding max diagonal element
	double result = fabs(diagonal_matrix[0][0]);
	for (int i = 1; i < size; ++i) {
		if (result < fabs(diagonal_matrix[i][i])) {
			result = fabs(diagonal_matrix[i][i]);
		}
	}

	return sqrt(result);
}

Matrix Matrix::get_diagonal_matrix() const {
	int im = 0;
	int jm = 0;
	// Copying matrix
	Matrix new_matrix = (*this);

	Matrix tmp_matrix(size);

	while (new_matrix.max_not_diagonal(im, jm) > eps) {
		double alpha = atan(2 * new_matrix[im][jm] / (new_matrix[im][im] - new_matrix[jm][jm])) / 2;

		Matrix matrix_T(size);
		for (int i = 0; i < matrix_T.size; ++i) {
			matrix_T[i][i] = 1;
		}

		matrix_T[im][im] = cos(alpha);
		matrix_T[jm][jm] = cos(alpha);
		matrix_T[im][jm] = -sin(alpha);
		matrix_T[jm][im] = sin(alpha);

		// T^T * A
		for (int i = 0; i < size; ++i) {
			for (int k = 0; k < size; ++k) {
				tmp_matrix[i][k] = 0;
				for (int j = 0; j < size; ++j) {
					tmp_matrix[i][k] += matrix_T[j][i] * new_matrix[j][k];
				}
			}
		}

		new_matrix = tmp_matrix * matrix_T;
	}

	return new_matrix;
}

double Matrix::max_not_diagonal(int& i, int& j) const {
	double res = -1;
	for (int k = 0; k < size; ++k) {
		for (int p = 0; p < size; ++p) {
			if (k == p)
				continue;

			if (res < fabs((*this)[k][p])) {
				res = fabs((*this)[k][p]);
				i = k;
				j = p;
			}
		}
	}

	return res;
}

/**
 * Simple iteration method for solving system of linear equations
 */
Vector Matrix::simple_iteration_method(const Vector& b, std::ostream& ostr) const {
	ostr << "                        | Discrepancy |  Error    |" << std::endl;
	ostr << "Itr |    tau   |  q     | Norm        |  Estimate |     x[1]      x[2]      x[3]      x[4]" << std::endl;
	Vector x = b;
	Vector prev_x(size);
	Vector next_x(size);
	int itr = 1;
	const double tau = 1.8 / get_third_norm();
	double discrepancy_norm;
	do {
		// find next x
		for (int i = 0; i < size; ++i) {
			double sum = 0;
			for (int j = 0; j < size; ++j) {
				sum += (*this)[i][j] * x[j];
			}

			next_x[i] = x[i] + tau * (b[i] - sum);
		}

		discrepancy_norm = (b - (*this) * next_x).find_third_norm();
		double q = find_transition_matrix_norm(itr, prev_x, x, next_x);
		double e = error_estimate(q, prev_x, x);

		prev_x = x;
		x = next_x;
		ostr << std::fixed << std::setprecision(PRECISION) << itr << " | " << tau << " | " << q << " | " <<
			discrepancy_norm << " | " << e << " | " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
		itr++;
	} while (discrepancy_norm > eps);

	return x;
}