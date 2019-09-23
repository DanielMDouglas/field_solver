#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <functional>

int wrap(int, bool, int);

// unit vectors
const std::vector <double> xhat = {1., 0., 0.};
const std::vector <double> yhat = {0., 1., 0.};
const std::vector <double> zhat = {0., 0., 1.};
const std::vector <double> zero_vector = {0., 0., 0.};

// handy vector initializer
std::vector <double> linspace(double, double, int);

// functionals for initializing boundary conditions 
std::function <double (double, double, double)> linear(double, double, double, double);
std::function <double (double, double, double)> constant(double);
std::function <double (double, double, double)> gaussian(double, double, double, double, double, double);

// functions of functionals
std::function <double (double, double, double)> func_sum(std::vector<std::function <double (double, double, double)>>);
std::function <double (double, double, double)> func_mul(std::vector<std::function <double (double, double, double)>>);
std::function <double (double, double, double)> func_rect(std::function <double (double, double, double)>);

// vector operators
std::vector <double> operator+(std::vector <double>, std::vector <double>);
std::vector <double> operator-(std::vector <double>, std::vector <double>);
std::vector <double> operator*(double, std::vector <double>);
std::vector <double> operator*(std::vector <double>, double);

// vector functions
double dot(std::vector <double>, std::vector <double>);
double mag(std::vector <double>);
std::vector <double> norm(std::vector <double>);

#endif
