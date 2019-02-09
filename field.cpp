#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "field.h"

field::field(std::vector<double> x, std::vector<double> y, double init_val = 0) {
  xSize = x.size();
  ySize = y.size();

  x_space = x;
  y_space = y;

  values = std::vector< double > (xSize*ySize, init_val);
}

field::field(std::vector<double> x, std::vector<double> y, double (*f)(double, double)) {
  xSize = x.size();
  ySize = y.size();

  x_space = x;
  y_space = y;

  values = std::vector< double > (xSize*ySize);

  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      set(i, j, (*f)(x[i], y[j]));
    }
  }
}

void field::set(int i, int j, double value) {
  values[ySize*i + j] = value;
}

double field::get(int i, int j) {
  return values[ySize*i + j];
}

void field::print_to_file(std::string filename) {
  std::ofstream outFile (filename.c_str());
  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      outFile << x_space[i] << ','
	      << y_space[j] << ','
	      << get(i, j) << '\n';
    }
  }
  outFile.close();
}

double squared_diff(field a, field b) {
  // return the squared differences of each field

  // first, check that they have the same shape
  try {
    if ( ( not a.xSize == b.xSize ) or ( not a.ySize == b.ySize ) ) {
      throw 20;
    }
    else {
      double sum = 0;
      int n_points = 0;
      for ( int i = 0; i < a.xSize; i++ ) {
	for ( int j = 0; j < a.ySize; j++ ) {
	  sum += pow(a.get(i,j) - b.get(i,j), 2);
	  n_points++;
	}
      }
      return sum/n_points;
    }
  }
  catch ( int e ) {
    std::cout << "ERROR: tried to compute squared difference of fields of different size!" << std::endl;
  }
  return 0;
}

double squared_diff(field a, field b, bool (*exclude)(int, int)) {
  // return the squared differences of each field

  // first, check that they have the same shape
  try {
    if ( ( not a.xSize == b.xSize ) or ( not a.ySize == b.ySize ) ) {
      throw 20;
    }
    else {
      double sum = 0;
      int n_points = 0;
      for ( int i = 0; i < a.xSize; i++ ) {
	for ( int j = 0; j < a.ySize; j++ ) {
	  if ( not exclude(i, j) ) {
	    sum += pow(a.get(i,j) - b.get(i,j), 2);
	    n_points++;
	  }
	}
      }
      return sum/n_points;
    }
  }
  catch ( int e ) {
    std::cout << "ERROR: tried to compute squared difference of fields of different size!" << std::endl;
  }
  return 0;
}
