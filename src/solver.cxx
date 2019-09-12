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
std::string startingQ = "none";
std::string outFileName = "final.dat";
std::string geom = "bulkPix";
int nIter = 100000;
int iterFreq = 100;
double tolerance = 1.e-15;
double w = 1;
int N = 0xdeadbeef;
double spacing = 0.01; // cm  
bool verbose = false;

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
    if ( optValue.str() == "-q" ) {
      argValue >> startingQ;
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
    if ( optValue.str() == "-f" ) {
      argValue >> iterFreq;
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
    if (optValue.str() == "-v" ) {
      verbose = true;
      opt--;
    }
  }

  std::cout << "####################################" << '\n'
	    << "# Using arguments: \n"
	    << "# inFile:              " << startingSol << '\n'
	    << "# outFile:             " << outFileName << '\n'
	    << "# geometry:            " << geom << '\n'
	    << "# max # iter:          " << nIter << '\n'
	    << "# reporting frequency: " << iterFreq << '\n'
	    << "# threshold:           " << tolerance << '\n'
            << "# N vertices:          " << N << '\n'
            << "# relaxation factor:   " << w << '\n'
	    << "# spacing:             " << spacing << '\n'
	    << "# verbose:             " << verbose << '\n'
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

  std::cout << "# "
	    << "Npoints X: " << nPointsX << '\t'
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
  std::vector <double> x_axis_staggered = linspace((b -> Xmin) - spacing/2,
						   (b -> Xmax) + spacing/2,
						   nPointsX + 1);
  std::vector <double> y_axis_staggered = linspace((b -> Ymin) - spacing/2,
						   (b -> Ymax) + spacing/2,
						   nPointsY + 1);
  std::vector <double> z_axis_staggered = linspace((b -> Zmin) - spacing/2,
						   (b -> Zmax) + spacing/2,
						   nPointsZ + 1);  

  // boundary voltage
  std::function <double (double, double, double)> bval_func;
  bval_func = std::function<double (double, double, double)> (std::bind(&boundary::boundary_value,
									b,
									std::placeholders::_1,
									std::placeholders::_2,
									std::placeholders::_3));
  bval = new field <double> (x_axis, y_axis, z_axis, bval_func);
			     

  // permittivity is defined everywhere
  // this gets permittivity from each volume, defaulting to 1 (vacuum)
  std::function <double (double, double, double)> epVal_func;
  epVal_func = std::function<double (double, double, double)> (std::bind(&boundary::permittivity,
									 b,
									 std::placeholders::_1,
									 std::placeholders::_2,
									 std::placeholders::_3));
  epVal = new field <double> (x_axis_staggered,
			      y_axis_staggered,
			      z_axis_staggered,
			      epVal_func);

  std::function <double (double, double, double)> sigVal_func;
  sigVal_func = std::function<double (double, double, double)> (std::bind(&boundary::conductivity,
									  b,
									  std::placeholders::_1,
									  std::placeholders::_2,
									  std::placeholders::_3));
  sigVal = new field <double> (x_axis_staggered,
			       y_axis_staggered,
			       z_axis_staggered,
			       sigVal_func);

  // is this point inside of the volume?
  std::function <bool (double, double, double)> is_in_vol_func;
  is_in_vol_func = std::function<bool (double, double, double)> (std::bind(&boundary::is_in_volume,
									   b,
									   std::placeholders::_1,
									   std::placeholders::_2,
									   std::placeholders::_3));
  is_in_volume = new field <bool> (x_axis, y_axis, z_axis, is_in_vol_func);

  // is this point inside of a conductor volume?
  std::function <bool (double, double, double)> is_dirichlet_func;
  is_dirichlet_func = std::function<bool (double, double, double)> (std::bind(&boundary::is_in_conductor,
									      b,
									      std::placeholders::_1,
									      std::placeholders::_2,
									      std::placeholders::_3));
  is_dirichlet = new field <bool> (x_axis, y_axis, z_axis, is_dirichlet_func);

  // is this point inside of a von neumann BC?
  std::function <bool (double, double, double)> is_vn_func;
  is_vn_func = std::function<bool (double, double, double)> (std::bind(&boundary::is_in_VN,
								       b,
								       std::placeholders::_1,
								       std::placeholders::_2,
								       std::placeholders::_3));
  is_von_neumann = new field <bool> (x_axis, y_axis, z_axis, is_vn_func);

  std::function <double (double, double, double)> vn_E_func;
  vn_E_func = std::function<double (double, double, double)> (std::bind(&boundary::Efield,
									b,
									std::placeholders::_1,
									std::placeholders::_2,
									std::placeholders::_3));
  von_neumann_dE = new field <double> (x_axis, y_axis, z_axis, vn_E_func);
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	von_neumann_dE -> set(i, j, k, (von_neumann_dE -> get(i, j, k))*spacing);
      }
    }
  }

  // charge distribution
  if ( startingQ == "none" ) {
    // start with 0 charge
    Q = new field <double> (x_axis, y_axis, z_axis, constant(0.));
  }
  else {
    Q = new field <double> (startingQ);
  }
  dQdt = new field <double> (x_axis, y_axis, z_axis, constant(0.));
  
  // field which will hold the iterative solution
  if ( startingSol == "none" ) {
    // if no starting solution is specified,
    // start by guessing a linear field
    double intercept = b -> boundary_value(0, 0, 0);
    double xSlope = (b -> boundary_value(b -> Xmax, 0, 0) - b -> boundary_value(b -> Xmin, 0, 0))/(b -> Xmax - b -> Xmin);
    double ySlope = (b -> boundary_value(0, b -> Ymax, 0) - b -> boundary_value(0, b -> Ymin, 0))/(b -> Ymax - b -> Ymin);
    double zSlope = (b -> boundary_value(0, 0, b -> Zmax) - b -> boundary_value(0, 0, b -> Zmin))/(b -> Zmax - b -> Zmin);
    // potential = new field <double> (x_axis, y_axis, z_axis,
    // 				    linear(intercept, xSlope, ySlope, zSlope));
    potential = new field <double> (x_axis, y_axis, z_axis,
				    constant(1));
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
	if ( is_dirichlet -> get(i, j, k) ) {
	  // set initial boundary values
	  potential -> set(i, j, k, bound -> boundary_value(x_axis[i],
							    y_axis[j],
							    z_axis[k]));
	}
	else {	    
	  // set permittivity coefficients
	  double e111 = epVal -> get(wrap(i+1, bound -> periodicX, nPointsX + 1),
				     wrap(j+1, bound -> periodicY, nPointsY + 1),
				     wrap(k+1, bound -> periodicZ, nPointsZ + 1));
	  double e011 = epVal -> get(wrap(i, bound -> periodicX, nPointsX + 1),
				     wrap(j+1, bound -> periodicY, nPointsY + 1),
				     wrap(k+1, bound -> periodicZ, nPointsZ + 1));
	  double e101 = epVal -> get(wrap(i+1, bound -> periodicX, nPointsX + 1),
				     wrap(j, bound -> periodicY, nPointsY + 1),
				     wrap(k+1, bound -> periodicZ, nPointsZ + 1));
	  double e110 = epVal -> get(wrap(i+1, bound -> periodicX, nPointsX + 1),
				     wrap(j+1, bound -> periodicY, nPointsY + 1),
				     wrap(k, bound -> periodicZ, nPointsZ + 1));
	  double e001 = epVal -> get(wrap(i, bound -> periodicX, nPointsX + 1),
				     wrap(j, bound -> periodicY, nPointsY + 1),
				     wrap(k+1, bound -> periodicZ, nPointsZ + 1));
	  double e010 = epVal -> get(wrap(i, bound -> periodicX, nPointsX + 1),
				     wrap(j+1, bound -> periodicY, nPointsY + 1),
				     wrap(k, bound -> periodicZ, nPointsZ + 1));
	  double e100 = epVal -> get(wrap(i+1, bound -> periodicX, nPointsX + 1),
				     wrap(j, bound -> periodicY, nPointsY + 1),
				     wrap(k, bound -> periodicZ, nPointsZ + 1));
	  double e000 = epVal -> get(wrap(i, bound -> periodicX, nPointsX + 1),
				     wrap(j, bound -> periodicY, nPointsY + 1),
				     wrap(k, bound -> periodicZ, nPointsZ + 1));
	  
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
	  double s111 = sigVal -> get(wrap(i+1, bound -> periodicX, nPointsX),
				      wrap(j+1, bound -> periodicY, nPointsY),
				      wrap(k+1, bound -> periodicZ, nPointsZ));
	  double s011 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j+1, bound -> periodicY, nPointsY),
				      wrap(k+1, bound -> periodicZ, nPointsZ));
	  double s101 = sigVal -> get(wrap(i+1, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k+1, bound -> periodicZ, nPointsZ));
	  double s110 = sigVal -> get(wrap(i+1, bound -> periodicX, nPointsX),
				      wrap(j+1, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));
	  double s001 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k+1, bound -> periodicZ, nPointsZ));
	  double s010 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j+1, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));
	  double s100 = sigVal -> get(wrap(i+1, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));
	  double s000 = sigVal -> get(wrap(i, bound -> periodicX, nPointsX),
				      wrap(j, bound -> periodicY, nPointsY),
				      wrap(k, bound -> periodicZ, nPointsZ));

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
  std::cout << "#\n"
	    << "# Iteration:"
	    << '\t'
	    << "Stepwise Difference Norm:"
	    << std::endl;
  for ( int iter = 0; iter < nIter; iter++ ) {
    // set tempGrid to the average of the neighboring cells
    for ( int i = 0; i < nPointsX; i++ ) {
      for ( int j = 0; j < nPointsY; j++ ) {
	for ( int k = 0; k < nPointsZ; k++ ) {
	  if ( (not is_dirichlet -> get(i, j, k))
	       and (not is_von_neumann -> get(i, j, k)) ) {
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

	    sum -= von_neumann_dE -> get(wrap(i-1, bound -> periodicX, nPointsX),
					 wrap(j, bound -> periodicY, nPointsY),
					 wrap(k, bound -> periodicZ, nPointsZ));
	    sum -= von_neumann_dE -> get(wrap(i+1, bound -> periodicX, nPointsX),
					 wrap(j, bound -> periodicY, nPointsY),
					 wrap(k, bound -> periodicZ, nPointsZ));
	    sum -= von_neumann_dE -> get(wrap(i, bound -> periodicX, nPointsX),
					 wrap(j-1, bound -> periodicY, nPointsY),
					 wrap(k, bound -> periodicZ, nPointsZ));
	    sum -= von_neumann_dE -> get(wrap(i, bound -> periodicX, nPointsX),
					 wrap(j+1, bound -> periodicY, nPointsY),
					 wrap(k, bound -> periodicZ, nPointsZ));
	    sum -= von_neumann_dE -> get(wrap(i, bound -> periodicX, nPointsX),
					 wrap(j, bound -> periodicY, nPointsY),
					 wrap(k-1, bound -> periodicZ, nPointsZ));
	    sum -= von_neumann_dE -> get(wrap(i, bound -> periodicX, nPointsX),
					 wrap(j, bound -> periodicY, nPointsY),
					 wrap(k+1, bound -> periodicZ, nPointsZ));

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

    if ( iter % iterFreq == 0 ) {
      // print out the squared_diff between current iteration and previous iteration
      double sqsm = squared_sum(resid, is_dirichlet);
      if ( not verbose ) {
	std::cout << '\r';
      }
      std::cout << iter << '\t' << '\t'
		<< sqsm << "        " << std::flush;
		// << std::endl;
      if ( verbose ) {
	std::cout << '\n';
      }

      if ( sqsm < tolerance ) {
	std::cout << "#\n"
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

  // set the value of the potential where it is von Neumann
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( is_von_neumann -> get(i, j, k) ) {
	  potential -> set(i, j, k, potential -> get(i, j, k-1) - von_neumann_dE -> get(i, j, k));
	  
	  // // each von Neumann grid point should have only one
	  // // non-von Neumann neighbor, so go through each of them...
	  // if ( (not is_von_neumann -> get(wrap(i+1, bound -> periodicX, nPointsX), j, k))
	  //      and (is_in_volume -> get(wrap(i+1, bound -> periodicX, nPointsX), j, k)) ) {
	  //   // std::cout << "VN point " << i << " " << j << " " << k << '\n'
	  //   // 	      << i+1 << " " << j << " " << k << std::endl;
	  //   potential -> set(i, j, k, von_neumann_dE -> get(i, j, k)
	  // 		     - potential -> get(wrap(i+1, bound -> periodicX, nPointsX), j, k));
	  // }
	  // else if ( (not is_von_neumann -> get(wrap(i-1, bound -> periodicX, nPointsX), j, k))
	  // 	    and (is_in_volume -> get(wrap(i-1, bound -> periodicX, nPointsX), j, k)) ) {
	  //   // std::cout << "VN point " << i << " " << j << " " << k << '\n'
	  //   // 	      << i-1 << " " << j << " " << k << std::endl;
	  //   potential -> set(i, j, k, von_neumann_dE -> get(i, j, k)
	  // 		     - potential -> get(wrap(i-1, bound -> periodicX, nPointsX), j, k));
	  // }
	  // else if ( (not is_von_neumann -> get(i, wrap(j+1, bound -> periodicY, nPointsY), k))
	  // 	    and (is_in_volume -> get(i, wrap(j+1, bound -> periodicY, nPointsY), k)) ) {
	  //   // std::cout << "VN point " << i << " " << j << " " << k << '\n'
	  //   // 	      << i << " " << j+1 << " " << k << std::endl;
	  //   potential -> set(i, j, k, von_neumann_dE -> get(i, j, k)
	  // 		     - potential -> get(i, wrap(j+1, bound -> periodicY, nPointsY), k));
	  // }
	  // else if ( (not is_von_neumann -> get(i, wrap(j-1, bound -> periodicY, nPointsY), k))
	  // 	    and (is_in_volume -> get(i, wrap(j-1, bound -> periodicY, nPointsY), k)) ) {
	  //   // std::cout << "VN point " << i << " " << j << " " << k << '\n'
	  //   // 	      << i << " " << j-1 << " " << k << std::endl;
	  //   potential -> set(i, j, k, von_neumann_dE -> get(i, j, k)
	  // 		     - potential -> get(i, wrap(j-1, bound -> periodicY, nPointsY), k));
	  // }
	  // else if ( (not is_von_neumann -> get(i, j, wrap(k+1, bound -> periodicZ, nPointsZ)))
	  // 	    and (is_in_volume -> get(i, j, wrap(k+1, bound -> periodicZ, nPointsZ))) ) {
	  //   // std::cout << "VN point " << i << " " << j << " " << k << '\n'
	  //   // 	      << i << " " << j << " " << k+1 << std::endl;
	  //   potential -> set(i, j, k, von_neumann_dE -> get(i, j, k)
	  // 		     - potential -> get(i, j, wrap(k+1, bound -> periodicZ, nPointsZ)));
	  // }
	  // else if ( (not is_von_neumann -> get(i, j, wrap(k-1, bound -> periodicZ, nPointsZ)))
	  // 	    and (is_in_volume -> get(i, j, wrap(k-1, bound -> periodicZ, nPointsZ))) ) {
	  //   // std::cout << "VN point " << i << " " << j << " " << k << '\n'
	  //   // 	      << i << " " << j << " " << k-1 << std::endl;
	  //   potential -> set(i, j, k, von_neumann_dE -> get(i, j, k)
	  // 		     - potential -> get(i, j, wrap(k-1, bound -> periodicZ, nPointsZ)));
	  // }
	}
      }
    }
  }

  std::cout << "# final stepwise difference: "
	    << squared_sum(resid, is_dirichlet)
	    << std::endl
	    << "#\n";
}

void solver::accum_charge()
{
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( not ( is_dirichlet -> get(i, j, k) or is_von_neumann -> get(i, j, k) ) ) {
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
  double sqsum = squared_sum(dQdt, is_dirichlet);
  std::cout << sqsum << std::endl;
  std::cout << squared_sum(Q, is_dirichlet) << std::endl; 
  int t = 0;
  double current_threshold = 1.e-32;
  
  while ( squared_sum(dQdt, is_dirichlet) > current_threshold ) {

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
		   + (dQdt -> get(i, j, k))*dt*1.e2);
	}
      }
    }
    sqsum = squared_sum(dQdt, is_dirichlet);
    solve_static();
    accum_charge();

    if ( sig_int ) {
      // check for ^C
      break;
    }
    std::cout << "Squared sum of dQ/dt " << sqsum << std::endl;
    std::cout << "Squared sum of Q " << squared_sum(Q, is_dirichlet) << std::endl;

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
  thisSolver.epVal -> print_to_file("perm_map.dat");
  thisSolver.potential -> print_to_file("initial.dat");
  
  // thisSolver.solve_charge();
  thisSolver.solve_static();
  thisSolver.potential -> print_to_file(outFileName);
  thisSolver.Q -> print_to_file("Q.dat");
    
  return 0;
}
