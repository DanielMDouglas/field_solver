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

solver::solver(boundary * b, int N, double spacing)
{
  bound = b;

  // set up grid parameters
  if ( N == 0xdeadbeef ) {
    nPointsX = (b -> Xmax - b -> Xmin)/spacing;
    nPointsY = (b -> Ymax - b -> Ymin)/spacing;
    nPointsZ = (b -> Zmax - b -> Zmin)/spacing;
  }
  else {
    nPointsX = N;
    nPointsY = N;
    nPointsZ = N;
    
    // hopefully these are the same for x, y, and z
    spacing = (b -> Xmax - b -> Xmin)/nPointsX;
  }
  if ( nPointsX == 0 ) {
    nPointsX++;
  }
  if ( nPointsY == 0 ) {
    nPointsY++;
  }
  if ( nPointsZ == 0 ) {
    nPointsZ++;
  }

  std::cout << "Npoints X: " << nPointsX << '\t'
	    << "Npoints Y: " << nPointsY << '\t'
	    << "Npoints Z: " << nPointsZ << '\t'
	    << std::endl;
  // over-relaxation: w > 1
  // does not work yet!
  // double t = cos(pi/nPointsX) + cos(pi/nPointsY) + cos(pi/nPointsZ);
  // double w = (8 - sqrt(64 - 16*t*t))/(t*t);

  std::vector <double> x_axis = linspace(b -> Xmin, b -> Xmax, nPointsX);
  std::vector <double> y_axis = linspace(b -> Ymin, b -> Ymax, nPointsY);
  std::vector <double> z_axis = linspace(b -> Zmin, b -> Zmax, nPointsZ);

  // staggered axes are for the permittivity grid
  // I think there might need to be another point (nPoints) for periodic boundary purposes
  std::vector <double> x_axis_staggered = linspace((b -> Xmin) + spacing/2, (b -> Xmax) + spacing/2, nPointsX);
  std::vector <double> y_axis_staggered = linspace((b -> Ymin) + spacing/2, (b -> Ymax) + spacing/2, nPointsY);
  std::vector <double> z_axis_staggered = linspace((b -> Zmin) + spacing/2, (b -> Zmax) + spacing/2, nPointsZ);  

  // boundary voltage
  bval = new field <double> (b, x_axis, y_axis, z_axis, "val");
  // permittivity is defined everywhere
  // this gets permittivity from each volume, defaulting to 1
  epVal = new field <double> (b, x_axis_staggered, y_axis_staggered, z_axis_staggered, "permittivity");
  sigVal = new field <double> (b, x_axis_staggered, y_axis_staggered, z_axis_staggered, "conductivity");
  // is this point inside of a conductor volume?
  is_b = new field <bool> (b, x_axis, y_axis, z_axis, "bool");
  // charge distribution
  // start with 0 charge
  Q = new field <double> (x_axis, y_axis, z_axis, constant(0.));
  dQdt = new field <double> (x_axis, y_axis, z_axis, constant(0.));

  // field which will hold the iterative solution
  if ( startingSol == "none" ) {
    // if no starting solution is specified,
    // start by guessing a linear field
    double intercept = b -> boundary_value(0, 0, 0);
    double xSlope = (b -> boundary_value(b -> Xmax, 0, 0) - b -> boundary_value(b -> Xmin, 0, 0))/(b -> Xmax - b -> Xmin);
    double ySlope = (b -> boundary_value(0, b -> Ymax, 0) - b -> boundary_value(0, b -> Ymin, 0))/(b -> Ymax - b -> Ymin);
    double zSlope = (b -> boundary_value(0, 0, b -> Zmax) - b -> boundary_value(0, 0, b -> Zmin))/(b -> Zmax - b -> Zmin);
    potential = new field <double> (x_axis, y_axis, z_axis,
				    linear(intercept, xSlope, ySlope, zSlope));
  }
  else {
    // otherwise, starting potential is loaded from a file
    potential = new field <double> (startingSol);
  }

  // initialize temporary fields
  resid = new field <double> (x_axis, y_axis, z_axis, constant(0.));
  a0 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  a1 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  a2 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  a3 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  a4 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  a5 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  a6 = new field <double> (x_axis, y_axis, z_axis, constant(1.));

  b0 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  b1 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  b2 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  b3 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  b4 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  b5 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
  b6 = new field <double> (x_axis, y_axis, z_axis, constant(1.));
    
  // initialize things that don't change over iterations  
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( is_b -> get(i, j, k) ) {
	  // set initial boundary values
	  potential -> set(i, j, k, bound -> boundary_value(x_axis[i],
						 y_axis[j],
						 z_axis[k]));
	}
	else {	    
	  // set permittivity coefficients
	  double e111 = epVal -> get(wrap(i, bound -> periodicX, nPointsX),
				     wrap(j, bound -> periodicY, nPointsY),
				     wrap(k, bound -> periodicZ, nPointsZ));
	  double e011 = epVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				     wrap(j, bound -> periodicY, nPointsY),
				     wrap(k, bound -> periodicZ, nPointsZ));
	  double e101 = epVal -> get(wrap(i, bound -> periodicX, nPointsX),
				     wrap(j-1, bound -> periodicY, nPointsY),
				     wrap(k, bound -> periodicZ, nPointsZ));
	  double e110 = epVal -> get(wrap(i, bound -> periodicX, nPointsX),
				     wrap(j, bound -> periodicY, nPointsY),
				     wrap(k-1, bound -> periodicZ, nPointsZ));
	  double e001 = epVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				     wrap(j-1, bound -> periodicY, nPointsY),
				     wrap(k, bound -> periodicZ, nPointsZ));
	  double e010 = epVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				     wrap(j, bound -> periodicY, nPointsY),
				     wrap(k-1, bound -> periodicZ, nPointsZ));
	  double e100 = epVal -> get(wrap(i, bound -> periodicX, nPointsX),
				     wrap(j-1, bound -> periodicY, nPointsY),
				     wrap(k-1, bound -> periodicZ, nPointsZ));
	  double e000 = epVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				     wrap(j-1, bound -> periodicY, nPointsY),
				     wrap(k-1, bound -> periodicZ, nPointsZ));
	  
	  // corresponds to central point
	  a0 -> set(i, j, k, 0.75*(e111 + e011 + e101 + e110
				   + e001 + e010 + e100 + e000));
	  // corresponds to +1 in k direction
	  a1 -> set(i, j, k, 0.25*(e111 + e011 + e101 + e001));
	  // corresponds to +1 in j direction
	  a2 -> set(i, j, k, 0.25*(e111 + e011 + e110 + e010));
	  // corresponds to +1 in i direction
	  a3 -> set(i, j, k, 0.25*(e111 + e101 + e110 + e100));
	  // corresponds to -1 in k direction
	  a4 -> set(i, j, k, 0.25*(e110 + e010 + e100 + e000));
	  // corresponds to -1 in j direction
	  a5 -> set(i, j, k, 0.25*(e101 + e001 + e100 + e000));
	  // corresponds to -1 in i direction
	  a6 -> set(i, j, k, 0.25*(e011 + e001 + e010 + e000));
	  	  
	  // set conductivity coefficients
	  double s111 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));
	  double s011 = sigVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));
	  double s101 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j-1, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));
	  double s110 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k-1, bound -> periodicZ, nPointsZ));
	  double s001 = sigVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				      wrap(j-1, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));
	  double s010 = sigVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k-1, bound -> periodicZ, nPointsZ));
	  double s100 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j-1, bound -> periodicY, nPointsY),
				      wrap(k-1, bound -> periodicZ, nPointsZ));
	  double s000 = sigVal -> get(wrap(i-1, bound -> periodicX, nPointsX),
				      wrap(j-1, bound -> periodicY, nPointsY),
				      wrap(k-1, bound -> periodicZ, nPointsZ));

	  // corresponds to central point
	  b0 -> set(i, j, k, 0.75*(s111 + s011 + s101 + s110 + s001 + s010 + s100 + s000));
	  // corresponds to +1 in k direction
	  b1 -> set(i, j, k, 0.25*(s111 + s011 + s101 + s001));
	  // corresponds to +1 in j direction
	  b2 -> set(i, j, k, 0.25*(s111 + s011 + s110 + s010));
	  // corresponds to +1 in i direction
	  b3 -> set(i, j, k, 0.25*(s111 + s101 + s110 + s100));
	  // corresponds to -1 in k direction
	  b4 -> set(i, j, k, 0.25*(s110 + s010 + s100 + s000));
	  // corresponds to -1 in j direction
	  b5 -> set(i, j, k, 0.25*(s101 + s001 + s100 + s000));
	  // corresponds to -1 in i direction
	  b6 -> set(i, j, k, 0.25*(s011 + s001 + s010 + s000));
	}
      }
    }
  }
}

void solver::solve_static()
{
  std::cout << "Iteration:"
	    << '\t'
	    << "Stepwise Difference Norm:"
	    << std::endl;
  for ( int iter = 0; iter < nIter; iter++ ) {
    // set tempGrid to the average of the neighboring cells
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  if ( not is_b -> get(i, j, k) ) {
	    double sum = 0;

	    sum += a1 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						       wrap(j, bound -> periodicY, nPointsY),
						       wrap(k+1, bound -> periodicZ, nPointsZ));
	    sum += a2 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						       wrap(j+1, bound -> periodicY, nPointsY),
						       wrap(k, bound -> periodicZ, nPointsZ));
	    sum += a3 -> get(i, j, k)*potential -> get(wrap(i+1, bound -> periodicX, nPointsX),
						       wrap(j, bound -> periodicY, nPointsY),
						       wrap(k, bound -> periodicZ, nPointsZ));
	    sum += a4 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						       wrap(j, bound -> periodicY, nPointsY),
						       wrap(k-1, bound -> periodicZ, nPointsZ));
	    sum += a5 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						       wrap(j-1, bound -> periodicY, nPointsY),
						       wrap(k, bound -> periodicZ, nPointsZ));
      	    sum += a6 -> get(i, j, k)*potential -> get(wrap(i-1, bound -> periodicX, nPointsX),
						       wrap(j, bound -> periodicY, nPointsY),
						       wrap(k, bound -> periodicZ, nPointsZ));

	    sum += (Q -> get(i, j, k))/ep0;
	    sum /= a0 -> get(i, j, k);
	    sum -= potential -> get(i, j, k);
	    resid -> set(i, j, k, sum);
	  }
	}
      }
    }

    if ( sig_int ) {
      // check for ^C
      std::cout << std::endl
		<< "# Recieved SIGTERM! Saving and quitting..."
		<< std::endl;
      break;
    }

    if ( iter % 100 == 0 ) {
      // print out the squared_diff between current iteration and previous iteration
      double sqsm = squared_sum(resid, is_b);
      std::cout << '\r'
		<< iter << '\t' << '\t'
		<< sqsm << "        " << std::flush;
		// << std::endl;
      if ( sqsm < tolerance ) {
	std::cout << std::endl
		  << "# potential converged after "
		  << iter
		  << " iterations"
		  << std::endl;
	break;
      }
    }
    
    for ( int i = 0; i < potential -> xSize; i++ ) {
      for ( int j = 0; j < potential -> ySize; j++ ) {
    	for ( int k = 0; k < potential -> zSize; k++ ) {
    	  double fac = w;
    	  potential -> set(i, j, k,
			   potential -> get(i, j, k) + fac*resid -> get(i, j, k));
    	}
      }
    }

  }
  
  std::cout << "# final stepwise difference: "
	    << squared_sum(resid, is_b)
	    << std::endl
	    << std::endl;
}

void solver::accum_charge()
{
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( not is_b -> get(i, j, k) ) {
	  double sum = 0;

	  sum += b1 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						     wrap(j, bound -> periodicY, nPointsY),
						     wrap(k+1, bound -> periodicZ, nPointsZ));
	  sum += b2 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						     wrap(j+1, bound -> periodicY, nPointsY),
						     wrap(k, bound -> periodicZ, nPointsZ));
	  sum += b3 -> get(i, j, k)*potential -> get(wrap(i+1, bound -> periodicX, nPointsX),
						     wrap(j, bound -> periodicY, nPointsY),
						     wrap(k, bound -> periodicZ, nPointsZ));
	  sum += b4 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						     wrap(j, bound -> periodicY, nPointsY),
						     wrap(k-1, bound -> periodicZ, nPointsZ));
	  sum += b5 -> get(i, j, k)*potential -> get(wrap(i, bound -> periodicX, nPointsX),
						     wrap(j-1, bound -> periodicY, nPointsY),
						     wrap(k, bound -> periodicZ, nPointsZ));
	  sum += b6 -> get(i, j, k)*potential -> get(wrap(i-1, bound -> periodicX, nPointsX),
						     wrap(j, bound -> periodicY, nPointsY),
						     wrap(k, bound -> periodicZ, nPointsZ));
	  sum -= b0 -> get(i, j, k)*potential -> get(i, j, k);

	  dQdt -> set(i, j, k, sum);
	}
      }
    }
  }
}

void solver::solve_charge()
{
  // Iteratively solve the static potential,
  // then evaluate div(J) to find charge accumulation,
  // then update charge distribution and solve new potential
  
  solve_static();
  accum_charge();
  double sqsum = squared_sum(dQdt, is_b);
  std::cout << sqsum << std::endl;
  std::cout << squared_sum(Q, is_b) << std::endl; 
  int t = 0;
  double current_threshold = 1.e-8;
  
  while ( squared_sum(dQdt, is_b) > current_threshold ) {

    std::stringstream VfileName, QfileName;
    VfileName << "V_series/" << t << ".dat";
    QfileName << "Q_series/" << t << ".dat";
    
    potential -> print_to_file(VfileName.str());
    Q -> print_to_file(QfileName.str());
    
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  Q -> set(i, j, k,
		   (Q -> get(i, j, k))
		   + (dQdt -> get(i, j, k))*dt*1.e-7);
	}
      }
    }
    sqsum = squared_sum(dQdt, is_b);
    solve_static();
    accum_charge();

    if ( sig_int ) {
      // check for ^C
      break;
    }
    std::cout << "Squared sum of dQ/dt " << sqsum << std::endl;
    std::cout << "Squared sum of Q " << squared_sum(Q, is_b) << std::endl;

    t++;
  }
}

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);

  // sig handling
  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = term;
  sigaction(SIGINT, &action, NULL);
  
  boundary * detector = new boundary (geom);
  solver thisSolver (detector, N, spacing);
  thisSolver.sigVal -> print_to_file("cond_map.dat");
  thisSolver.solve_charge();
  // thisSolver.solve_static();
  thisSolver.potential -> print_to_file(outFileName);
  thisSolver.Q -> print_to_file("Q.dat");
  
  // thisSolver.solve_static();
  // thisSolver.potential -> print_to_file(outFileName);
  // thisSolver.accum_charge();
  // thisSolver.dQdt -> print_to_file("drho.dat");
  
  return 0;
}
