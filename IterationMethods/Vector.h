#pragma once
#include <iostream>

class Vector {
	int size{};
	double* arr{};

public:
	Vector() = delete;
	explicit Vector(int n);
	Vector(const Vector& other);
	Vector& operator=(Vector other);
	Vector(Vector&& other) noexcept;
	~Vector();

	friend void swap(Vector& lhs, Vector& rhs) noexcept;

	double& operator[](int index) const;
	void input(std::istream& fin) const;
	void print(std::ostream& ostr) const;
	double find_first_norm() const;
	double find_second_norm() const;
	double find_third_norm() const;
	Vector operator-(const Vector& other) const;
	Vector operator+(const Vector& other) const;
	Vector operator*(double number) const;
	double operator*(const Vector& other) const;
};