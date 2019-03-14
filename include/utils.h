#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <functional>

// unit vectors
const std::vector <double> xhat = {1., 0., 0.};
const std::vector <double> yhat = {0., 1., 0.};
const std::vector <double> zhat = {0., 0., 1.};

// handy vector initializer
std::vector <double> linspace(double, double, int);

// functionals for initializing boundary conditions 
std::function<double (double, double, double)> linear(double, double, double, double);
std::function<double (double, double, double)> constant(double);

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
