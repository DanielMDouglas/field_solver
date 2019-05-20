#include <cassert>
#include <vector>
#include <math.h>
#include <iostream>

#include "utils.h"

int wrap(int index, bool periodicity, int maxIndex) {
  int newIndex;
  if ( index == -1 ) {
    if ( periodicity ) {
      newIndex = maxIndex - 1;
    }
    else {
      newIndex = 0xdeadbeef; // not nice, but magic number that scalarfield::get will handle
    }
  }
  else if ( index == maxIndex ) {
    if ( periodicity ) {
      newIndex = 0;
    }
    else {
      newIndex = 0xdeadbeef;
    }
  }
  else {
    newIndex = index;
  }

  return newIndex; 
}

std::vector <double> linspace(double start, double end, int nPoints)
{
  std::vector <double> result (nPoints);
  double delta = (end - start)/(nPoints - 1);
  for ( int i = 0; i < nPoints; i++ ) {
    result[i] = start + i*delta;
  }
  return result;
}

std::function <double (double, double, double)> linear(double intercept,
						      double xSlope,
						      double ySlope,
						      double zSlope)
{
  // intercept defines the value of the function at the origin (x = y = z = 0)
  // x/y/zSlope defines the change in the value w.r.t. x/y/z
  return [intercept, xSlope, ySlope, zSlope] (double x, double y, double z)
    {
      return intercept + xSlope*x + ySlope*y + zSlope*z;
    };
}

std::function <double (double, double, double)> constant(double C)
{
  return [C](double, double, double) { return C; };
}

std::function <double (double, double, double)>  gaussian(double x0,
							  double y0,
							  double z0,
							  double width,
							  double baseline,
							  double magnitude)
{
  return [x0, y0, z0, width, baseline, magnitude] (double x, double y, double z)
 	 {
	   return baseline + magnitude*exp(-1*(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))/(2*pow(width, 2)));
	 };
}

std::function <double (double, double, double)> func_sum(std::vector<std::function <double (double, double, double)>> funcList)
{
  return [funcList] (double x, double y, double z)
	 {
	   double sum = 0;
	   for ( std::function <double (double, double, double)> thisFunc : funcList ) {
	     sum += thisFunc(x, y, z);
	   }
	   return sum;
	 };
}

std::function <double (double, double, double)> func_mul(std::vector<std::function <double (double, double, double)>> funcList)
{
  return [funcList] (double x, double y, double z)
	 {
	   double prod = 1;
	   for ( std::function <double (double, double, double)> thisFunc : funcList ) {
	     prod *= thisFunc(x, y, z);
	   }
	   return prod;
	 };
}

std::function <double (double, double, double)> func_rect(std::function <double (double, double, double)> func)
{
  // return a "rectified" version of the argument
  // that is, the same function with negative values mapped to 0
  return [func] (double x, double y, double z)
	 {
	   double result = func(x, y, z);
	   if ( result < 0 ) {
	     result = 0;
	   }
	   return result;
	 };
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
