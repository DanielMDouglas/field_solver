#include <cassert>
#include <vector>
#include <math.h>

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

std::vector <double> operator+(std::vector <double> first, std::vector <double> second)
{
  assert ( first.size() == second.size() );

  std::vector <double> result;
  result.reserve(first.size());

  for ( uint i = 0; i < first.size(); i++ ) {
    result.push_back(first[i] + second[i]);
  }

  return result;
}

std::vector <double> operator-(std::vector <double> first, std::vector <double> second)
{
  assert ( first.size() == second.size() );

  std::vector <double> result;
  result.reserve(first.size());

  for ( uint i = 0; i < first.size(); i++ ) {
    result.push_back(first[i] - second[i]);
  }

  return result;
}

std::vector <double> operator*(double scalar, std::vector <double> vect)
{
  std::vector <double> result;
  result.reserve(vect.size());

  for ( uint i = 0; i < vect.size(); i++ ) {
    result.push_back(scalar*vect[i]);
  }

  return result;
}

std::vector <double> operator*(std::vector <double> vect, double scalar)
{
  return scalar*vect;
}

double dot(std::vector <double> first, std::vector <double> second)
{
  assert ( first.size() == second.size() );
  
  double sum = 0;

  for ( uint i = 0; i < first.size(); i++ ) {
    sum += first[i] * second[i];
  }

  return sum;
}

double mag(std::vector <double> vect)
{
  return sqrt(dot(vect, vect));
}

std::vector <double> norm(std::vector <double> vect)
{
  return (1./mag(vect))*vect;
}
