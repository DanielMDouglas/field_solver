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
int N = 280;

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
  }

  std::cout << "Using arguments: \n"
	    << "inFile:            " << startingSol << '\n'
	    << "outFile:           " << outFileName << '\n'
	    << "geometry:          " << geom << '\n'
	    << "max # iter:        " << nIter << '\n'
	    << "threshold:         " << tolerance << '\n'
            << "N vertices:        " << N << '\n'
            << "relaxation factor: " << w << std::endl;
}

void term(int signum)
{
  sig_int = 1;
}

void solve_field(boundary b)
{
  // make some initial guess for the solution
  // assume a linear field
  
  int nPointsX = N;
  int nPointsY = N;
  int nPointsZ = N;  

  // double t = cos(pi/nPointsX) + cos(pi/nPointsY) + cos(pi/nPointsZ);
  // double w = (8 - sqrt(64 - 16*t*t))/(t*t);
  
  scalarField <double> bval (b, nPointsX, nPointsY, nPointsZ, "val");
  scalarField <bool> is_b (b, nPointsX, nPointsY, nPointsZ, "bool");

  std::cout << b.Xmin << '\t' << b.Xmax << '\n'
	    << b.Ymin << '\t' << b.Ymax << '\n'
	    << b.Zmin << '\t' << b.Zmax << std::endl;
  std::vector <double> x_axis = linspace(b.Xmin, b.Xmax, nPointsX);
  std::vector <double> y_axis = linspace(b.Ymin, b.Ymax, nPointsY);
  std::vector <double> z_axis = linspace(b.Zmin, b.Zmax, nPointsZ);

  scalarField <double> analytic_solution (x_axis, y_axis, z_axis, linear(0, 0, 0, 1000));
  analytic_solution.print_to_file("analytic.dat");
  
  scalarField <double> solution (x_axis, y_axis, z_axis, constant(0.));
  if ( startingSol == "none" ) {
    double intercept = b.boundary_value(0, 0, 0);
    double xSlope = (b.boundary_value(b.Xmax, 0, 0) - b.boundary_value(b.Xmin, 0, 0))/(b.Xmax - b.Xmin);
    double ySlope = (b.boundary_value(0, b.Ymax, 0) - b.boundary_value(0, b.Ymin, 0))/(b.Ymax - b.Ymin);
    double zSlope = (b.boundary_value(0, 0, b.Zmax) - b.boundary_value(0, 0, b.Zmin))/(b.Zmax - b.Zmin);

    std::cout << intercept << '\t' << xSlope << '\t' << ySlope << '\t' << zSlope << std::endl;
    intercept = 0;
    
    // solution = scalarField <double> (x_axis, y_axis, z_axis,
    // 				     linear(intercept, xSlope, ySlope, zSlope));
  }
  else {
    solution = scalarField <double> (startingSol);
  }

  scalarField <double> resid (x_axis, y_axis, z_axis, constant(0.));
  scalarField <double> weight (x_axis, y_axis, z_axis, constant(1.));
  
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
	  weight.set(i, j, k, 0.);
	  if ( not ( i == 0 ) ) {  
	    weight.set(i, j, k, weight.get(i, j, k) + 1);
	  }
	  else {
	    if ( b.periodicX ) {
	      weight.set(i, j, k, weight.get(i, j, k) + 1);
	    }
	  }
	  if ( not ( i == nPointsX-1 ) ) {
	    weight.set(i, j, k, weight.get(i, j, k) + 1);
	  }
	  else {
	    if ( b.periodicX ) {
	      weight.set(i, j, k, weight.get(i, j, k) + 1);
	    }
	  }
	  if ( not ( j == 0 ) ) {
	    weight.set(i, j, k, weight.get(i, j, k) + 1);
	  }
	  else {
	    if ( b.periodicY ) {
	      weight.set(i, j, k, weight.get(i, j, k) + 1);
	    }
	  }
	  if ( not ( j == nPointsY-1 ) ) {
	    weight.set(i, j, k, weight.get(i, j, k) + 1);
	  }
	  else {
	    if ( b.periodicY ) {
	      weight.set(i, j, k, weight.get(i, j, k) + 1);
	    }
	  }
	  if ( not ( k == 0 ) ) {
	    weight.set(i, j, k, weight.get(i, j, k) + 1);
	  }
	  else {
	    if ( b.periodicZ ) {
	      weight.set(i, j, k, weight.get(i, j, k) + 1);
	    }
	  }
	  if ( not ( k == nPointsZ-1 ) ) {
	    weight.set(i, j, k, weight.get(i, j, k) + 1);
	  }
	  else {
	    if ( b.periodicZ ) {
	      weight.set(i, j, k, weight.get(i, j, k) + 1);
	    }
	  }
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
		sum += solution.get(nPointsX-1, j, k);
	      }
	    }
	    else {
	      sum += solution.get(i-1, j, k);
	    }
	    if ( i == nPointsX-1 ) {
	      if ( b.periodicX ) {
		sum += solution.get(0, j, k);
	      }
	    }
	    else {
	      sum += solution.get(i+1, j, k);
	    }

	    // y boundary conditions
	    if ( j == 0 ) {
	      if ( b.periodicY ) {
		sum += solution.get(i, nPointsY-1, k);
	      }
	    }
	    else {
	      sum += solution.get(i, j-1, k);
	    }
	    if ( j == nPointsY-1 ) {
	      if ( b.periodicY ) {
		sum += solution.get(i, 0, k);
	      }
	    }
	    else {
	      sum += solution.get(i, j+1, k);
	    }

	    // z boundary conditions
	    if ( k == 0 ) {
	      if ( b.periodicZ ) {
		sum += solution.get(i, j, nPointsZ-1);
	      }
	    }
	    else {
	      sum += solution.get(i, j, k-1);
	    }
	    if ( k == nPointsZ-1 ) {
	      if ( b.periodicZ ) {
		sum += solution.get(i, j, 0);
	      }
	    }
	    else {
	      sum += solution.get(i, j, k+1);
	    }

	    sum /= weight.get(i, j, k);
	    sum -= solution.get(i, j, k);
	    resid.set(i, j, k, sum);
	  }
	}
      }
    }

    if ( sig_int ) {
      std::cout << std::endl << "Recieved SIGTERM! Saving and quitting..." << std::endl;
      break;
    }

    if ( iter % 10 == 0 ) {
      // print out the squared_diff between current iteration and previous iteration
      std::cout << iter << '\t'
		<< squared_sum(resid, is_b) << '\t'
		<< "boop" << '\t'
		<< std::endl;
      if ( squared_sum(resid, is_b) < tolerance ) {
	std::cout << "solution converged after " << iter << " iterations" << std::endl;
	break;
      }
    }
    
    // set solution to the tempGrid
    for ( int i = 0; i < solution.xSize*solution.ySize*solution.zSize; i++ ) {
      solution.values[i] = solution.values[i] + w*resid.values[i];
    }
  }
  
  std::cout << "final stepwise difference: " << squared_sum(resid, is_b) << std::endl;
  
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
