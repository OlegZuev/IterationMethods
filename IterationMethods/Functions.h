#pragma once
#include <iostream>

class Vector;
void print_header(std::ostream& ostr);

double find_transition_matrix_norm(const Vector& prev_x, const Vector& x, const Vector& next_x);

double error_estimate(double q, const Vector& x, const Vector& next_x);