#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <math.h>

#include "field.h"

double squared_sum(field <double> * a)
{
  // return the squared sum of a field

  double sum = 0;
  int n_points = 0;
  for ( int i = 0; i < a -> xSize; i++ ) {
    for ( int j = 0; j < a -> ySize; j++ ) {
      for ( int k = 0; k < a -> zSize; k++ ) {
	sum += pow(a -> get(i, j, k), 2);
	n_points++;
      }
    }
  }
  return sum/n_points;
}

double squared_sum(field <double> * a, field <bool> * exclude)
{
  // return the squared sum of a field where a boolean field is true

  double sum = 0;
  int n_points = 0;
  for ( int i = 0; i < a -> xSize; i++ ) {
    for ( int j = 0; j < a -> ySize; j++ ) {
      for ( int k = 0; k < a -> zSize; k++ ) {
	if ( not exclude -> get(i, j, k) ) {
	  sum += pow(a -> get(i, j, k), 2);
	  n_points++;
	}
      }
    }
  }
  return sum/n_points;
}

double squared_diff(field <double> * a, field <double> * b)
{
  // return the sum of squared differences of each field

  field <double> * temp = new field <double>(a -> x_space, a -> y_space, a -> z_space, 0.);
  for ( int i = 0; i < a -> xSize; i++ ) {
    for ( int j = 0; j < a -> ySize; j++ ) {
      for ( int k = 0; k < a -> zSize; k++ ) {
	temp -> set(i, j, k, (a -> get(i, j, k)) - (b -> get(i, j, k)));
      }
    }
  }
  return squared_sum(temp); 
}

double squared_diff(field <double> * a, field <double> * b, field <bool> * exclude)
{
  // return the sum of squared differences of each field where a boolean field is true

  field <double> * temp = new field <double> (a -> x_space, a -> y_space, a -> z_space, constant(0.));
  for ( int i = 0; i < (a -> xSize)*(a -> ySize)*(a -> zSize); i++ ) {
    temp -> values[i] = (a -> values[i]) - (b -> values[i]);	
  }
  return squared_sum(temp, exclude);
}
