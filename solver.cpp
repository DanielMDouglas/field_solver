#include <iostream>
#include <fstream>
#include <vector>

#include "field.h"

double xL = 1000;
double yL = 200;

// the size of each step (in x and y)
double ds = 2;
const int nPointsX = xL/ds;
const int nPointsY = yL/ds;

// number of iterations
// with ds = 2, 15000 seems to be a good number
int nIter = 5000;

// constant by which each timestep is multiplied
// should be 1./4 if dx == dy
double k = 1./4;

double Vmax = 100;

std::vector<double> linspace(double start, double end, int nPoints) {
  std::vector<double> result (nPoints);
  double delta = (end - start)/(nPoints - 1);
  for ( int i = 0; i < nPoints; i++ ) {
    result[i] = start + i*delta;
  }
  
  return result;
}

std::vector<double> x_axis = linspace(0., xL, nPointsX+1);
std::vector<double> y_axis = linspace(0., yL, nPointsY);

bool is_in_boundary(int i, int j) {
  if ( i == 0 ) {
    return true;
  }
  else if ( i == nPointsX ) {
    return true;
  }
  else if ( j == 0 ) {
    return true;
  }
  else if ( j == nPointsY-1) {
    return true;
  }
  else return false;
}

double boundary_value(int i, int j) {
  if ( i == 0 ) {
    return y_axis[j]/yL*Vmax;
  }
  else if ( i == nPointsX ) {
    return y_axis[j]/yL*Vmax;
  }
  else if ( j == 0 ) {
    return 0;
  }
  else if ( j == nPointsY-1) {
    return Vmax;
  }
  else return 0;
}

field V (x_axis, y_axis, Vmax/2);
field tempGrid (x_axis, y_axis, 0.);

double linear(double x, double y) {
  return y*Vmax/yL;
}

field analytic_solution(x_axis, y_axis, linear);  

int main() {

  // set initial boundary values
  // these values don't get updated
  for ( int i = 0; i < nPointsX+1; i++ ) {
    for ( int j = 0; j < nPointsY; j ++ ) {
      if ( is_in_boundary(i, j) ) {
	V.set(i, j, boundary_value(i, j));
      }
    }
  }

  for ( int iter = 0; iter < nIter; iter++ ) {
    // set tempGrid to the average 
    for ( int i = 0; i < nPointsX+1; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	if ( not is_in_boundary(i, j) ) {
	  tempGrid.set(i, j,
		       k*(V.get(i-1, j) +
			  V.get(i+1, j) +
			  V.get(i, j-1) +
			  V.get(i, j+1)));
	}
      }
    }
    // set V to the tempGrid
    for ( int i = 0; i < nPointsX+1; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	if ( not is_in_boundary(i, j) ) {
	  V.set(i, j,
		tempGrid.get(i,j));
	}
      }
    }
  }
  
  V.print_to_file("convergence.dat");
  analytic_solution.print_to_file("analytic.dat");
  
  return 0;
}
