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
std::string outFileName = "final.root";
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

void solve_field(boundary b)
{
  // set up grid parameters
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
    
  // over-relaxation: w > 1
  // does not work yet!
  // double t = cos(pi/nPointsX) + cos(pi/nPointsY) + cos(pi/nPointsZ);
  // double w = (8 - sqrt(64 - 16*t*t))/(t*t);
  
  std::vector <double> x_axis = linspace(b.Xmin, b.Xmax, nPointsX);
  std::vector <double> y_axis = linspace(b.Ymin, b.Ymax, nPointsY);
  std::vector <double> z_axis = linspace(b.Zmin, b.Zmax, nPointsZ);

  // staggered axes are for the permittivity grid
  // I think there might need to be another point (nPoints) for periodic boundary purposes
  std::vector <double> x_axis_staggered = linspace(b.Xmin + spacing/2, b.Xmax - spacing/2, nPointsX-1);
  std::vector <double> y_axis_staggered = linspace(b.Ymin + spacing/2, b.Ymax - spacing/2, nPointsY-1);
  std::vector <double> z_axis_staggered = linspace(b.Zmin + spacing/2, b.Zmax - spacing/2, nPointsZ-1);  

  // boundary voltage
  // field <double> bval (b, x_axis, y_axis, z_axis, "val");
  field <double> bval (b, x_axis, y_axis, z_axis, "val");
  // permittivity is defined everywhere
  // this gets permittivity from each volume, defaulting to 1
  field <double> epVal (b, x_axis_staggered, y_axis_staggered, z_axis_staggered, "permittivity");
  // is this point inside of a conductor volume?
  field <bool> is_b (b, x_axis, y_axis, z_axis, "bool");
  // charge distribution
  // start with 0 charge
  field <double> Q (x_axis, y_axis, z_axis, constant(0.));
  
  // field which will hold the iterative solution
  field <double> solution (x_axis, y_axis, z_axis, constant(0.));
  if ( startingSol == "none" ) {
    // if no starting solution is specified,
    // start by guessing a linear field
    double intercept = b.boundary_value(0, 0, 0);
    double xSlope = (b.boundary_value(b.Xmax, 0, 0) - b.boundary_value(b.Xmin, 0, 0))/(b.Xmax - b.Xmin);
    double ySlope = (b.boundary_value(0, b.Ymax, 0) - b.boundary_value(0, b.Ymin, 0))/(b.Ymax - b.Ymin);
    double zSlope = (b.boundary_value(0, 0, b.Zmax) - b.boundary_value(0, 0, b.Zmin))/(b.Zmax - b.Zmin);
    solution = field <double> (x_axis, y_axis, z_axis,
			       linear(intercept, xSlope, ySlope, zSlope));
  }
  else {
    // otherwise, starting solution is loaded from a file
    solution = field <double> (startingSol);
  }

  // initialize temporary fields
  field <double> resid (x_axis, y_axis, z_axis, constant(0.));
  field <double> weight (x_axis, y_axis, z_axis, constant(1.));
  field <double> a0 (x_axis, y_axis, z_axis, constant(1.));
  field <double> a1 (x_axis, y_axis, z_axis, constant(1.));
  field <double> a2 (x_axis, y_axis, z_axis, constant(1.));
  field <double> a3 (x_axis, y_axis, z_axis, constant(1.));
  field <double> a4 (x_axis, y_axis, z_axis, constant(1.));
  field <double> a5 (x_axis, y_axis, z_axis, constant(1.));
  field <double> a6 (x_axis, y_axis, z_axis, constant(1.));
  
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
	  // set permittivity coefficients
	  double e111 = epVal.get(wrap(i, b.periodicX, nPointsX),
	  			  wrap(j, b.periodicY, nPointsY),
	  			  wrap(k, b.periodicZ, nPointsZ));
	  double e011 = epVal.get(wrap(i-1, b.periodicX, nPointsX),
	  			  wrap(j, b.periodicY, nPointsY),
	  			  wrap(k, b.periodicZ, nPointsZ));
	  double e101 = epVal.get(wrap(i, b.periodicX, nPointsX),
	  			  wrap(j-1, b.periodicY, nPointsY),
	  			  wrap(k, b.periodicZ, nPointsZ));
	  double e110 = epVal.get(wrap(i, b.periodicX, nPointsX),
	  			  wrap(j, b.periodicY, nPointsY),
	  			  wrap(k-1, b.periodicZ, nPointsZ));
	  double e001 = epVal.get(wrap(i-1, b.periodicX, nPointsX),
	  			  wrap(j-1, b.periodicY, nPointsY),
	  			  wrap(k, b.periodicZ, nPointsZ));
	  double e010 = epVal.get(wrap(i-1, b.periodicX, nPointsX),
	  			  wrap(j, b.periodicY, nPointsY),
	  			  wrap(k-1, b.periodicZ, nPointsZ));
	  double e100 = epVal.get(wrap(i, b.periodicX, nPointsX),
	  			  wrap(j-1, b.periodicY, nPointsY),
	  			  wrap(k-1, b.periodicZ, nPointsZ));
	  double e000 = epVal.get(wrap(i-1, b.periodicX, nPointsX),
	  			  wrap(j-1, b.periodicY, nPointsY),
	  			  wrap(k-1, b.periodicZ, nPointsZ));
	  
	  // corresponds to central point
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

	    sum += a1.get(i, j, k)*solution.get(wrap(i, b.periodicX, nPointsX),
	    					wrap(j, b.periodicY, nPointsY),
	    					wrap(k+1, b.periodicZ, nPointsZ));
	    sum += a2.get(i, j, k)*solution.get(wrap(i, b.periodicX, nPointsX),
	    					wrap(j+1, b.periodicY, nPointsY),
	    					wrap(k, b.periodicZ, nPointsZ));
	    sum += a3.get(i, j, k)*solution.get(wrap(i+1, b.periodicX, nPointsX),
	    					wrap(j, b.periodicY, nPointsY),
	    					wrap(k, b.periodicZ, nPointsZ));
	    sum += a4.get(i, j, k)*solution.get(wrap(i, b.periodicX, nPointsX),
	    					wrap(j, b.periodicY, nPointsY),
	    					wrap(k-1, b.periodicZ, nPointsZ));
	    sum += a5.get(i, j, k)*solution.get(wrap(i, b.periodicX, nPointsX),
	    					wrap(j-1, b.periodicY, nPointsY),
	    					wrap(k, b.periodicZ, nPointsZ));
      	    sum += a6.get(i, j, k)*solution.get(wrap(i-1, b.periodicX, nPointsX),
	    					wrap(j, b.periodicY, nPointsY),
	    					wrap(k, b.periodicZ, nPointsZ));

	    sum += Q.get(i, j, k);
	    sum /= a0.get(i, j, k);
	    sum -= solution.get(i, j, k);
	    resid.set(i, j, k, sum);
	  }
	}
      }
    }

    if ( sig_int ) {
      // check for ^C
      std::cout << std::endl << "# Recieved SIGTERM! Saving and quitting..." << std::endl;
      break;
    }

    if ( iter % 50 == 0 ) {
      // print out the squared_diff between current iteration and previous iteration
      double sqsm = squared_sum(resid, is_b);
      std::cout << iter << '\t'
		<< sqsm << '\t'
		<< std::endl;
      // is this small enough?
      if ( sqsm < tolerance ) {
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
