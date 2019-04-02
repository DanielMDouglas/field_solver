#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <signal.h>
#include <cstring>

#include "solver.h"
#include "physics.h"

volatile int sig_int = 0;
std::string startingSol = "none";
std::string outFileName = "final.dat";
std::string geom = "bulkPix";
int nIter = 100000;
double tolerance = 1.e-15;
double w = 1;
int N = 0xdeadbeef; //previously 280
double spacing = 0.01; // cm  

void handleOpts(int argc, char const * argv[])
{
  int opt = 0;
  while ( opt < argc ) {
    std::stringstream optValue;
    std::stringstream argValue;
    optValue << argv[++opt];
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
    if ( optValue.str() == "-w" ) {
      argValue >> w;
    }
    if ( optValue.str() == "-N" ) {
      argValue >> N;
    }
    if ( optValue.str() == "-s" ) {
      argValue >> spacing;
    }
  }

  std::cout << "####################################" << '\n'
	    << "# Using arguments: \n"
	    << "# inFile:            " << startingSol << '\n'
	    << "# outFile:           " << outFileName << '\n'
	    << "# geometry:          " << geom << '\n'
	    << "# max # iter:        " << nIter << '\n'
	    << "# threshold:         " << tolerance << '\n'
            << "# N vertices:        " << N << '\n'
            << "# relaxation factor: " << w << '\n'
	    << "# spacing:           " << spacing << '\n'
	    << "####################################" << std::endl;
}

void term(int signum)
{
  sig_int = 1;
}

void solve_field(boundary b)
{
  // make some initial guess for the solution
  // assume a linear field

  int nPointsX;
  int nPointsY;
  int nPointsZ;
  if ( N == 0xdeadbeef ) {
    nPointsX = (b.Xmax - b.Xmin)/spacing;
    nPointsY = (b.Ymax - b.Ymin)/spacing;
    nPointsZ = (b.Zmax - b.Zmin)/spacing;
  }
  else {
    nPointsX = N;
    nPointsY = N;
    nPointsZ = N;
    
    // hopefully these are the same for x, y, and z
    double spacing = (b.Xmax - b.Xmin)/nPointsX;
  }
    
  // double t = cos(pi/nPointsX) + cos(pi/nPointsY) + cos(pi/nPointsZ);
  // double w = (8 - sqrt(64 - 16*t*t))/(t*t);
  
  std::vector <double> x_axis = linspace(b.Xmin, b.Xmax, nPointsX);
  std::vector <double> y_axis = linspace(b.Ymin, b.Ymax, nPointsY);
  std::vector <double> z_axis = linspace(b.Zmin, b.Zmax, nPointsZ);

  std::vector <double> x_axis_staggered = linspace(b.Xmin + spacing/2, b.Xmax - spacing/2, nPointsX-1);
  std::vector <double> y_axis_staggered = linspace(b.Ymin + spacing/2, b.Ymax - spacing/2, nPointsY-1);
  std::vector <double> z_axis_staggered = linspace(b.Zmin + spacing/2, b.Zmax - spacing/2, nPointsZ-1);  

  scalarField <double> bval (b, x_axis, y_axis, z_axis, "val");
  scalarField <double> epVal (b, x_axis_staggered, y_axis_staggered, z_axis_staggered, "permittivity");
  scalarField <bool> is_b (b, x_axis, y_axis, z_axis, "bool");
  
  scalarField <double> analytic_solution (x_axis, y_axis, z_axis, linear(0, 0, 0, 1));
  analytic_solution.print_to_file("analytic.dat");
  
  scalarField <double> solution (x_axis, y_axis, z_axis, constant(0.));
  if ( startingSol == "none" ) {
    double intercept = b.boundary_value(0, 0, 0);
    double xSlope = (b.boundary_value(b.Xmax, 0, 0) - b.boundary_value(b.Xmin, 0, 0))/(b.Xmax - b.Xmin);
    double ySlope = (b.boundary_value(0, b.Ymax, 0) - b.boundary_value(0, b.Ymin, 0))/(b.Ymax - b.Ymin);
    double zSlope = (b.boundary_value(0, 0, b.Zmax) - b.boundary_value(0, 0, b.Zmin))/(b.Zmax - b.Zmin);
    
    solution = scalarField <double> (x_axis, y_axis, z_axis,
    				     linear(intercept, xSlope, ySlope, zSlope));
  }
  else {
    solution = scalarField <double> (startingSol);
  }

  scalarField <double> resid (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> weight (x_axis, y_axis, z_axis, constant(1.));
  scalarField <double> a0 (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> a1 (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> a2 (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> a3 (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> a4 (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> a5 (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> a6 (x_axis, y_axis, z_axis, constant(0.));
  
  // initialize things that don't change over iterations
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( is_b.get(i, j, k) ) {
	  // set initial boundary values
	  solution.set(i, j, k, b.boundary_value(x_axis[i],
						 y_axis[j],
						 z_axis[k]));
	}
	else {
	  double e111 = epVal.get(i, j, k);
	  double e011 = epVal.get(i-1, j, k);
	  double e101 = epVal.get(i, j-1, k);
	  double e110 = epVal.get(i, j, k-1);
	  double e001 = epVal.get(i-1, j-1, k);
	  double e010 = epVal.get(i-1, j, k-1);
	  double e100 = epVal.get(i, j-1, k-1);
	  double e000 = epVal.get(i-1, j-1, k-1);
	  
	  a0.set(i, j, k, 0.75*(e111 + e011 + e101 + e110
				+ e001 + e010 + e100 + e000));
	  // corresponds to +1 in k direction
	  a1.set(i, j, k, 0.25*(e111 + e011 + e101 + e001));
	  // corresponds to +1 in j direction
	  a2.set(i, j, k, 0.25*(e111 + e011 + e110 + e010));
	  // corresponds to +1 in i direction
	  a3.set(i, j, k, 0.25*(e111 + e101 + e110 + e100));
	  // corresponds to -1 in k direction
	  a4.set(i, j, k, 0.25*(e110 + e010 + e100 + e000));
	  // corresponds to -1 in j direction
	  a5.set(i, j, k, 0.25*(e101 + e001 + e100 + e000));
	  // corresponds to -1 in i direction
	  a6.set(i, j, k, 0.25*(e011 + e001 + e010 + e000));
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
      
	    // x boundary conditions
	    if ( i == 0 ) {
	      if ( b.periodicX ) {
		sum += a6.get(i, j, k)*solution.get(nPointsX-1, j, k);
	      }
	    }
	    else {
	      sum += a6.get(i, j, k)*solution.get(i-1, j, k);
	    }
	    if ( i == nPointsX-1 ) {
	      if ( b.periodicX ) {
		sum += a3.get(i, j, k)*solution.get(0, j, k);
	      }
	    }
	    else {
	      sum += a3.get(i, j, k)*solution.get(i+1, j, k);
	    }

	    // y boundary conditions
	    if ( j == 0 ) {
	      if ( b.periodicY ) {
		sum += a5.get(i, j, k)*solution.get(i, nPointsY-1, k);
	      }
	    }
	    else {
	      sum += a5.get(i, j, k)*solution.get(i, j-1, k);
	    }
	    if ( j == nPointsY-1 ) {
	      if ( b.periodicY ) {
		sum += a2.get(i, j, k)*solution.get(i, 0, k);
	      }
	    }
	    else {
	      sum += a2.get(i, j, k)*solution.get(i, j+1, k);
	    }

	    // z boundary conditions
	    if ( k == 0 ) {
	      if ( b.periodicZ ) {
		sum += a4.get(i, j, k)*solution.get(i, j, nPointsZ-1);
	      }
	    }
	    else {
	      sum += a4.get(i, j, k)*solution.get(i, j, k-1);
	    }
	    if ( k == nPointsZ-1 ) {
	      if ( b.periodicZ ) {
		sum += a1.get(i, j, k)*solution.get(i, j, 0);
	      }
	    }
	    else {
	      sum += a1.get(i, j, k)*solution.get(i, j, k+1);
	    }

	    // sum /= weight.get(i, j, k);
	    // sum -= a0.get(i, j, k)*solution.get(i, j, k);
	    sum /= a0.get(i, j, k);
	    sum -= solution.get(i, j, k);
	    resid.set(i, j, k, sum);
	  }
	}
      }
    }

    if ( sig_int ) {
      std::cout << std::endl << "# Recieved SIGTERM! Saving and quitting..." << std::endl;
      break;
    }

    if ( iter % 10 == 0 ) {
      // print out the squared_diff between current iteration and previous iteration
      std::cout << iter << '\t'
		<< squared_sum(resid, is_b) << '\t'
		<< squared_diff(solution, analytic_solution, is_b)  << '\t'
		<< std::endl;
      if ( squared_sum(resid, is_b) < tolerance ) {
	std::cout << "# solution converged after " << iter << " iterations" << std::endl;
	break;
      }
    }
    
    // set solution to the tempGrid
    // for ( int i = 0; i < solution.xSize*solution.ySize*solution.zSize; i++ ) {
    //   solution.values[i] = solution.values[i] + w*resid.values[i];
    // }
    for ( int i = 0; i < solution.xSize; i++ ) {
      for ( int j = 0; j < solution.ySize; j++ ) {
    	for ( int k = 0; k < solution.zSize; k++ ) {
    	  double fac = 1;
    	  if ( z_axis[k] > 1.0 ) {
    	    fac = w;
    	  }
    	  solution.set(i, j, k,
    		       solution.get(i, j, k) + fac*resid.get(i, j, k));
    	}
      }
    }

  }
  
  std::cout << "# final stepwise difference: " << squared_sum(resid, is_b) << std::endl;
  
  solution.print_to_file(outFileName);
}

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);

  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = term;
  sigaction(SIGINT, &action, NULL);
  
  boundary detector (geom);
  solve_field(detector);
  
  return 0;
}
