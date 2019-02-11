#include <iostream>
#include <fstream>
#include <vector>

#include "field.h"

double xL = 50;
double yL = 50;
double zL = 50;

// the size of each step (in x and y)
double ds = 1;
const int nPointsX = xL/ds;
const int nPointsY = yL/ds;
const int nPointsZ = zL/ds;

// number of iterations
// with ds = 2, 15000 seems to be a good number
int nIter = 10000;

// constant by which each timestep is multiplied
// should be 1./4 if dx == dy
// using 1./6 for 3 dimensions
double weight_f = 1./6;

double Vmax = 100;

std::vector<double> linspace(double start, double end, int nPoints) {
  std::vector<double> result (nPoints);
  double delta = (end - start)/(nPoints - 1);
  for ( int i = 0; i < nPoints; i++ ) {
    result[i] = start + i*delta;
  }
  
  return result;
}

std::vector <double> x_axis = linspace(0., xL, nPointsX);
std::vector <double> y_axis = linspace(0., yL, nPointsY);
std::vector <double> z_axis = linspace(0., zL, nPointsZ);

bool is_in_boundary(int i, int j, int k) {
  if ( i == 0 ) {
    return true;
  }
  else if ( i == nPointsX-1 ) {
    return true;
  }
  else if ( j == 0 ) {
    return true;
  }
  else if ( j == nPointsY-1 ) {
    return true;
  }
  else if ( k == 0 ) {
    return true;
  }
  else if ( k == nPointsZ-1 ) {
    return true;
  }
  else return false;
}

double boundary_value(int i, int j, int k) {
  if ( i == 0 ) {
    return y_axis[j]/yL*Vmax;
  }
  else if ( i == nPointsX-1 ) {
    return y_axis[j]/yL*Vmax;
  }
  if ( k == 0 ) {
    return y_axis[j]/yL*Vmax;
  }
  else if ( k == nPointsZ-1 ) {
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

field V (x_axis, y_axis, z_axis, Vmax/2);
field tempGrid (x_axis, y_axis, z_axis, 0.);

double linear(double x, double y, double z) {
  // the analytic solution
  // used to initialize the analytic field
  return y*Vmax/yL;
}

field analytic_solution(x_axis, y_axis, z_axis, linear);  

int main() {

  // set initial boundary values
  // these values don't get updated
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) { 
	if ( is_in_boundary(i, j, k) ) {
	  V.set(i, j, k, boundary_value(i, j, k));
	}
      }
    }
  }

  for ( int iter = 0; iter < nIter; iter++ ) {
    // set tempGrid to the average 
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  if ( not is_in_boundary(i, j, k) ) {
	    tempGrid.set(i, j, k,
			 weight_f*(V.get(i-1, j, k) +
				   V.get(i+1, j, k) +
				   V.get(i, j-1, k) +
				   V.get(i, j+1, k) +
				   V.get(i, j, k-1) +
				   V.get(i, j, k+1)));
	  }
	}
      }
    }

    if ( iter % 100 == 0 ) {
      // print out the squared_diff between current iteration
      // and the analytic solution
      std::cout << iter << '\t'
		<< squared_diff(V, analytic_solution, is_in_boundary) << '\t'
		<< squared_diff(V, tempGrid, is_in_boundary) << '\t'
		<< std::endl;
      if ( squared_diff(V, tempGrid, is_in_boundary) < 1.e-7 ) {
	break;
      }
    }
      
    // set V to the tempGrid
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  if ( not is_in_boundary(i, j, k) ) {
	    V.set(i, j, k,
		  tempGrid.get(i, j, k));
	  }
	}
      }
    }
  }

  std::cout << "final stepwise difference: " << squared_diff(V, tempGrid, is_in_boundary) << std::endl;
  
  V.print_to_file("convergence.dat");
  analytic_solution.print_to_file("analytic.dat");
  
  return 0;
}
