#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <signal.h>
#include <cstring>

#include "solver.h"

volatile int sig_int = 0;
std::string startingSol = "none";
std::string outFileName = "final.dat";
std::string geom = "bulkPix";
int nIter = 100000;
double tolerance = 1.e-15;

void handleOpts(int argc, char const *argv[])
{
  int opt = 1;
  std::stringstream optValue;
  std::stringstream argValue;
  while ( opt < argc ) {
    optValue << argv[opt];
    argValue << argv[++opt];
    if ( optValue.str() == "-i" ) {
      argValue >> startingSol;
    }
    if ( optValue.str() == "-o" ) {
      argValue >> outFileName;
    }
    if ( optValue.str() == "-g" ) {
      argValue >> geom;
    }
    if ( optValue.str() == "-n" ) {
      argValue >> nIter;
    }
    if ( optValue.str() == "-t" ) {
      argValue >> tolerance;
    }
    opt++;
  }

  std::cout << "Using arguments: \n"
	    << "inFile:     " << startingSol << '\n'
	    << "outFile:    " << outFileName << '\n'
	    << "geometry:   " << geom << '\n'
	    << "max # iter: " << nIter << '\n'
	    << "threshold:  " << tolerance << std::endl;
}

void term(int signum)
{
  sig_int = 1;
}

void solve_field(boundary b) {
  // make some initial guess for the solution
  // assume a linear field
  
  int nPointsX = 200;
  int nPointsY = 200;
  int nPointsZ = 200;  
  
  scalarField bval = scalarField(b, nPointsX, nPointsY, nPointsZ, "val");
  scalarField is_b = scalarField(b, nPointsX, nPointsY, nPointsZ, "bool");
  
  std::vector <double> x_axis = linspace(b.Xmin, b.Xmax, nPointsX);
  std::vector <double> y_axis = linspace(b.Ymin, b.Ymax, nPointsY);
  std::vector <double> z_axis = linspace(b.Zmin, b.Zmax, nPointsZ);

  scalarField solution = scalarField(x_axis, y_axis, z_axis, constant(0.));
  if ( startingSol == "none" ) {
    double intercept = b.boundary_value(0, 0, 0);
    double xSlope = (b.boundary_value(b.Xmax, 0, 0) - b.boundary_value(b.Xmin, 0, 0))/(b.Xmax - b.Xmin);
    double ySlope = (b.boundary_value(0, b.Ymax, 0) - b.boundary_value(0, b.Ymin, 0))/(b.Ymax - b.Ymin);
    double zSlope = (b.boundary_value(0, 0, b.Zmax) - b.boundary_value(0, 0, b.Zmin))/(b.Zmax - b.Zmin);

    solution = scalarField(x_axis, y_axis, z_axis,
			   linear(intercept, xSlope, ySlope, zSlope));
  }
  else {
    solution = scalarField(startingSol);
  }

  scalarField tempGrid = scalarField(x_axis, y_axis, z_axis, constant(0.));
  
  // set initial boundary values
  // these values don't get updated
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( is_b.get(i, j, k) ) {
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
	  if ( not is_b.get(i, j, k) ) {
	    double sum = 0;
	    double weight = 0;
	    if ( not ( i == 0 ) ) {  
	      sum += solution.get(i-1, j, k);
	      weight += 1;
	    }
	    else {
	      if ( b.periodicX ) {
		sum += solution.get(nPointsX-1, j, k);
		weight += 1;
	      }
	    }
	    if ( not ( i == nPointsX-1 ) ) {
	      sum += solution.get(i+1, j, k);
	      weight += 1;
	    }
	    else {
	      if ( b.periodicX ) {
		sum += solution.get(0, j, k);
		weight += 1;
	      }
	    }
	    if ( not ( j == 0 ) ) {
	      sum += solution.get(i, j-1, k);
	      weight += 1;
	    }
	    else {
	      if ( b.periodicY ) {
		sum += solution.get(i, nPointsY-1, k);
		weight += 1;
	      }
	    }
	    if ( not ( j == nPointsY-1 ) ) {
	      sum += solution.get(i, j+1, k);
	      weight += 1;
	    }
	    else {
	      if ( b.periodicY ) {
		sum += solution.get(i, 0, k);
		weight += 1;
	      }
	    }
	    if ( not ( k == 0 ) ) {
	      sum += solution.get(i, j, k-1);
	      weight += 1;
	    }
	    else {
	      if ( b.periodicZ ) {
		sum += solution.get(i, j, nPointsZ-1);
		weight += 1;
	      }
	    }
	    if ( not ( k == nPointsZ-1 ) ) {
	      sum += solution.get(i, j, k+1);
	      weight += 1;
	    }
	    else {
	      if ( b.periodicZ ) {
		sum += solution.get(i, j, 0);
		weight += 1;
	      }
	    }
	    sum /= weight;
	    tempGrid.set(i, j, k, sum);
	  }
	}
      }
    }

    if ( sig_int ) {
      std::cout << std::endl << "Recieved SIGTERM! Saving and quitting..." << std::endl;
      break;
    }
    
    if ( iter % 100 == 0 ) {
      // print out the squared_diff between current iteration and previous iteration
      std::cout << iter << '\t'
		<< squared_diff(solution, tempGrid, is_b) << '\t'
		<< std::endl;
      if ( squared_diff(solution, tempGrid, is_b) < tolerance ) {
	std::cout << "solution converged after " << iter << " iterations" << std::endl;
	break;
      }
    }
      
    // set solution to the tempGrid
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  if ( not is_b.get(i, j, k) ) {
	    solution.set(i, j, k,
		  tempGrid.get(i, j, k));
	  }
	}
      }
    }
  }

  std::cout << "final stepwise difference: " << squared_diff(solution, tempGrid, is_b) << std::endl;
  
  solution.print_to_file(outFileName);
}

int main(int argc, char const *argv[]) {
  handleOpts(argc, argv);

  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = term;
  sigaction(SIGINT, &action, NULL);
  
  boundary detector (geom);
  solve_field(detector);
  
  return 0;
}
