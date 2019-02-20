#include "utils.h"

std::vector <double> linspace(double start, double end, int nPoints)
{
  std::vector <double> result (nPoints);
  double delta = (end - start)/(nPoints - 1);
  for ( int i = 0; i < nPoints; i++ ) {
    result[i] = start + i*delta;
  }

  return result;
}

std::function<double (double, double, double)> linear(double intercept,
						      double xSlope,
						      double ySlope,
						      double zSlope)
{
  // intercept defines the value of the function at the origin
  // x/y/zSlope defines the change in the value w.r.t. x/y/z
  return [intercept, xSlope, ySlope, zSlope] (double x, double y, double z)
    {
      return intercept + xSlope*x + ySlope*y + zSlope*z;
    };
}

std::function<double (double, double, double)> constant(double C)
{
  return [C](double, double, double) { return C; };
}
