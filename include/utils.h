#include <vector>
#include <functional>

std::vector <double> linspace(double, double, int);

std::function<double (double, double, double)> linear(double, double, double, double);

std::function<double (double, double, double)> constant(double);
