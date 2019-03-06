#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <functional>

const std::vector <double> xhat = {1., 0., 0.};
const std::vector <double> yhat = {0., 1., 0.};
const std::vector <double> zhat = {0., 0., 1.};

std::vector <double> linspace(double, double, int);

std::function<double (double, double, double)> linear(double, double, double, double);

std::function<double (double, double, double)> constant(double);

std::vector <double> operator+(std::vector <double>, std::vector <double>);
std::vector <double> operator-(std::vector <double>, std::vector <double>);
std::vector <double> operator*(double, std::vector <double>);
std::vector <double> operator*(std::vector <double>, double);

double dot(std::vector <double>, std::vector <double>);
double mag(std::vector <double>);
std::vector <double> norm(std::vector <double>);

#endif
