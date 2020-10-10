#include "Vector.h"
#include <utility>
#include <cmath>

/**
 * Constructor with allocating memory for Vector
 */
Vector::Vector(int n) {
	size = n;
	arr = new double[size] {0};
}

/**
 * Copy constructor
 */
Vector::Vector(const Vector& other) : Vector(other.size) {
	for (int i = 0; i < size; ++i) {
		arr[i] = other.arr[i];
	}
}

/**
 * Copy assignment operator
 */
Vector& Vector::operator=(Vector other) {
	swap(*this, other);

	return *this;
}

/**
 * Move constructor
 */
Vector::Vector(Vector&& other) noexcept {
	swap(*this, other);
}

/**
 * Destructor
 */
Vector::~Vector() {
	delete[] arr;
}

/**
 * Swap
 */
void swap(Vector& lhs, Vector& rhs) noexcept {
	using std::swap;

	swap(lhs.size, rhs.size);
	swap(lhs.arr, rhs.arr);
}

/**
 * Indexer
 */
double& Vector::operator[](const int index) const {
	return arr[index];
}

/**
 * Input Vector from "fin" stream
 */
void Vector::input(std::istream& fin) const {
	for (int i = 0; i < size; ++i) {
		fin >> arr[i];
	}
}

/**
 * Print Vector into "ostr" stream
 */
void Vector::print(std::ostream& ostr) const {
	for (int i = 0; i < size; ++i) {
		ostr << arr[i] << " ";
	}

	ostr << std::endl;
}

/**
 * Find the third norm (Euclid) of vector
 */
double Vector::find_third_norm() const {
	double res = 0;
	for (int i = 0; i < size; i++) {
		res += (*this)[i] * (*this)[i];
	}

	return sqrt(res);
}

/**
 * subtraction of vectors
 */
Vector Vector::operator-(const Vector& other) const {
	Vector new_vector(size);
	for (int i = 0; i < size; ++i) {
		new_vector[i] = (*this)[i] - other[i];
	}

	return new_vector;
}

/**
 * sum of vectors
 */
Vector Vector::operator+(const Vector& other) const {
	Vector new_vector(size);
	for (int i = 0; i < size; ++i) {
		new_vector[i] = (*this)[i] + other[i];
	}

	return new_vector;
}

Vector Vector::operator*(const double number) const {
	Vector new_vector(size);
	for (int i = 0; i < size; i++) {
		new_vector[i] = arr[i] * number;
	}

	return new_vector;
}
