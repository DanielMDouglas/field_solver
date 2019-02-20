#include <iostream>
#include <fstream>
#include <vector>

#include "solver.h"

void solve_field(boundary b) {
  // make some initial guess for the solution
  // assume a linear field
  int nPointsX = 100;
  int nPointsY = 100;
  int nPointsZ = 100;

  double weight_f = 1./6;
  int nIter = 10000;
  
  field bval = field(b, nPointsX, nPointsY, nPointsZ, "val");
  field is_b = field(b, nPointsX, nPointsY, nPointsZ, "bool");
  
  double intercept = b.boundary_value(0, 0, 0);
  double xSlope = (b.boundary_value(b.Xmax, 0, 0) - b.boundary_value(b.Xmin, 0, 0))/(b.Xmax - b.Xmin);
  double ySlope = (b.boundary_value(0, b.Ymax, 0) - b.boundary_value(0, b.Ymin, 0))/(b.Ymax - b.Ymin);
  double zSlope = (b.boundary_value(0, 0, b.Zmax) - b.boundary_value(0, 0, b.Zmin))/(b.Zmax - b.Zmin);

  std::vector <double> x_axis = linspace(b.Xmin, b.Xmax, nPointsX);
  std::vector <double> y_axis = linspace(b.Ymin, b.Ymax, nPointsY);
  std::vector <double> z_axis = linspace(b.Zmin, b.Zmax, nPointsZ);

  field solution = field(x_axis, y_axis, z_axis, linear(intercept, xSlope, ySlope, zSlope));
  field tempGrid = field(x_axis, y_axis, z_axis, constant(0.));
  
  solution.print_to_file("initial.dat");

  // set initial boundary values
  // these values don't get updated
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) { 
	if ( b.is_in_boundary(x_axis[i],
			      y_axis[j],
			      z_axis[k]) ) {
	  solution.set(i, j, k, b.boundary_value(x_axis[i],
						 y_axis[j],
						 z_axis[k]));
	}
      }
    }
  }

  for ( int iter = 0; iter < nIter; iter++ ) {
    // set tempGrid to the average of the neighboring cells
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  if ( not b.is_in_boundary(x_axis[i],
				    y_axis[j],
				    z_axis[k]) ) {
	    tempGrid.set(i, j, k,
			 weight_f*(solution.get(i-1, j, k) +
				   solution.get(i+1, j, k) +
				   solution.get(i, j-1, k) +
				   solution.get(i, j+1, k) +
				   solution.get(i, j, k-1) +
				   solution.get(i, j, k+1)));
	  }
	}
      }
    }

    if ( iter % 100 == 0 ) {
      // print out the squared_diff between current iteration and previous iteration
      std::cout << iter << '\t'
		<< squared_diff(solution, tempGrid, is_b) << '\t'
		<< std::endl;
      if ( squared_diff(solution, tempGrid, is_b) < 1.e-7 ) {
	std::cout << "solution converged after " << iter << " iterations" << std::endl;
	break;
      }
    }
      
    // set solution to the tempGrid
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  if ( not b.is_in_boundary(x_axis[i], y_axis[j], z_axis[k]) ) {
	    solution.set(i, j, k,
		  tempGrid.get(i, j, k));
	  }
	}
      }
    }
  }

  std::cout << "final stepwise difference: " << squared_diff(solution, tempGrid, is_b) << std::endl;
  
  solution.print_to_file("final.dat");
}

int main() {

  boundary detector;
  solve_field(detector);
  
  return 0;
}
