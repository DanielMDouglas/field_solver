#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <math.h>

#include "scalar_field.h"

double squared_sum(scalarField <double> a)
{
  // return the squared sum of a field

  double sum = 0;
  int n_points = 0;
  for ( int i = 0; i < a.xSize; i++ ) {
    for ( int j = 0; j < a.ySize; j++ ) {
      for ( int k = 0; k < a.zSize; k++ ) {
	sum += pow(a.get(i, j, k), 2);
	n_points++;
      }
    }
  }
  return sum/n_points;
}

double squared_sum(scalarField <double> a, scalarField <bool> exclude)
{
  // return the squared sum of a field

  double sum = 0;
  int n_points = 0;
  for ( int i = 0; i < a.xSize; i++ ) {
    for ( int j = 0; j < a.ySize; j++ ) {
      for ( int k = 0; k < a.zSize; k++ ) {
	if ( not exclude.get(i, j, k) ) {
	  sum += pow(a.get(i, j, k), 2);
	  n_points++;
	}
      }
    }
  }
  return sum/n_points;
}

double squared_diff(scalarField <double> a, scalarField <double> b)
{
  // return the squared differences of each field

  scalarField <double> temp = scalarField <double>(a.x_space, a.y_space, a.z_space, 0.);
  for ( int i = 0; i < a.xSize; i++ ) {
    for ( int j = 0; j < a.ySize; j++ ) {
      for ( int k = 0; k < a.zSize; k++ ) {
	temp.set(i, j, k, a.get(i, j, k) - b.get(i, j, k));
      }
    }
  }
  return squared_sum(temp); 
}

double squared_diff(scalarField <double> a, scalarField <double> b, scalarField <bool> exclude)
{
  // return the squared differences of each field

  scalarField <double> temp (a.x_space, a.y_space, a.z_space, constant(0.));
  for ( int i = 0; i < a.xSize*a.ySize*a.zSize; i++ ) {
    temp.values[i] = a.values[i] - b.values[i];	
  }
  return squared_sum(temp, exclude);
}
