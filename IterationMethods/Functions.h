#pragma once

class Vector;
double find_third_norm_vector(double* x, int n);

double find_transition_matrix_norm(int itr, const Vector& prev_x, const Vector& x, const Vector& next_x);

double error_estimate(double q, const Vector& prev_x, const Vector& x);