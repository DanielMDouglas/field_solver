#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <signal.h>
#include <cstring>
#include <thread>

// #ifdef _OPENMP
#include <omp.h>
// #endif

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
double initSpacing = 0.01; // cm  
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
      argValue >> initSpacing;
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
	    << "# spacing:             " << initSpacing << '\n'
	    << "# verbose:             " << verbose << '\n'
	    << "####################################" << std::endl;
}

void term(int signum)
{
  sig_int = 1;
}

solver::solver(boundary * b, int N, double initSpacing)
{
  bound = b;
  spacing = initSpacing;

  // set up grid parameters
  if ( N == 0xdeadbeef ) {
    nPointsX = (bound -> Xmax - bound -> Xmin)/spacing + 1;
    nPointsY = (bound -> Ymax - bound -> Ymin)/spacing + 1;
    nPointsZ = (bound -> Zmax - bound -> Zmin)/spacing + 1;
  }
  else {
    nPointsX = N;
    nPointsY = N;
    nPointsZ = N;
    
    // hopefully these are the same for x, y, and z
    spacing = (bound -> Xmax - bound -> Xmin)/(nPointsX - 1);
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
	    << std::endl
	    << "# "
	    << "Spacing: " << spacing
	    << std::endl;
  // over-relaxation: w > 1
  // does not work yet!
  // double t = cos(pi/nPointsX) + cos(pi/nPointsY) + cos(pi/nPointsZ);
  // double w = (8 - sqrt(64 - 16*t*t))/(t*t);

  initialize_axes();

  initialize_geometry_fields();
  
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
    double intercept = bound -> boundary_value(0, 0, 0);
    double xSlope = (bound -> boundary_value(bound -> Xmax, 0, 0) - bound -> boundary_value(bound -> Xmin, 0, 0))/(bound -> Xmax - bound -> Xmin);
    double ySlope = (bound -> boundary_value(0, bound -> Ymax, 0) - bound -> boundary_value(0, bound -> Ymin, 0))/(bound -> Ymax - bound -> Ymin);
    double zSlope = (bound -> boundary_value(0, 0, bound -> Zmax) - bound -> boundary_value(0, 0, bound -> Zmin))/(bound -> Zmax - bound -> Zmin);
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
  initialize_permittivity_and_conductivity();
  
  // initialize things that don't change over iterations  
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( is_boundary -> get(i, j, k) ) {
	  // set initial boundary values
	  potential -> set(i, j, k, bound -> boundary_value(x_axis[i],
							    y_axis[j],
							    z_axis[k]));
	}
      }
    }
  }

  set_VN();
}

void solver::upscale(int scaleFactor)
{
  nPointsX = scaleFactor*(nPointsX - 1) + 1;
  nPointsY = scaleFactor*(nPointsY - 1) + 1;
  nPointsZ = scaleFactor*(nPointsZ - 1) + 1;
  
  spacing /= scaleFactor;

  std::cout << "# "
	    << "Npoints X: " << nPointsX << '\t'
	    << "Npoints Y: " << nPointsY << '\t'
	    << "Npoints Z: " << nPointsZ << '\t'
	    << std::endl
    	    << "# "
	    << "Spacing: " << spacing
	    << std::endl;

  initialize_axes();

  initialize_geometry_fields();

  Q = Q -> upscale(scaleFactor);

  potential = potential -> upscale(scaleFactor);

  // initialize temporary fields
  resid = new field <double> (x_axis, y_axis, z_axis, constant(0.));

  initialize_permittivity_and_conductivity();

  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( is_boundary -> get(i, j, k) ) {
	  // set initial boundary values
	  potential -> set(i, j, k, bound -> boundary_value(x_axis[i],
							    y_axis[j],
							    z_axis[k]));
	}
      }
    }
  }
  
  set_VN();
}

void solver::initialize_axes()
{
  // needs
  //     boundary * b
  //     int nPointsX
  //     int nPointsY
  //     int nPointsZ
  // to be set already
  
  x_axis = linspace(bound -> Xmin, bound -> Xmax, nPointsX);
  y_axis = linspace(bound -> Ymin, bound -> Ymax, nPointsY);
  z_axis = linspace(bound -> Zmin, bound -> Zmax, nPointsZ);

  x_axis_staggered = linspace((bound -> Xmin) - spacing/2,
			      (bound -> Xmax) + spacing/2,
			      nPointsX + 1);
  y_axis_staggered = linspace((bound -> Ymin) - spacing/2,
			      (bound -> Ymax) + spacing/2,
			      nPointsY + 1);
  z_axis_staggered = linspace((bound -> Zmin) - spacing/2,
			      (bound -> Zmax) + spacing/2,
			      nPointsZ + 1);
}

void solver::initialize_geometry_fields()
{
  // needs
  //     boundary * bound
  //     std::vector <double> x_axis
  //     std::vector <double> y_axis
  //     std::vector <double> z_axis
  //     std::vector <double> x_axis_staggered
  //     std::vector <double> y_axis_staggered
  //     std::vector <double> z_axis_staggered
  // to be set already
  
  // boundary voltage
  std::function <double (double, double, double)> bval_func;
  bval_func = std::function<double (double, double, double)> (std::bind(&boundary::boundary_value,
									bound,
									std::placeholders::_1,
									std::placeholders::_2,
									std::placeholders::_3));
  bval = new field <double> (x_axis, y_axis, z_axis, bval_func);

    // permittivity is defined everywhere
  // this gets permittivity from each volume, defaulting to 1 (vacuum)
  std::function <double (double, double, double)> epVal_func;
  epVal_func = std::function<double (double, double, double)> (std::bind(&boundary::permittivity,
									 bound,
									 std::placeholders::_1,
									 std::placeholders::_2,
									 std::placeholders::_3));
  epVal = new field <double> (x_axis_staggered,
			      y_axis_staggered,
			      z_axis_staggered,
			      epVal_func);

  std::function <double (double, double, double)> sigVal_func;
  sigVal_func = std::function<double (double, double, double)> (std::bind(&boundary::conductivity,
									  bound,
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
									   bound,
									   std::placeholders::_1,
									   std::placeholders::_2,
									   std::placeholders::_3));
  is_in_volume = new field <bool> (x_axis, y_axis, z_axis, is_in_vol_func);

  // is this point inside of a conductor volume?
  std::function <bool (double, double, double)> is_dirichlet_func;
  is_dirichlet_func = std::function<bool (double, double, double)> (std::bind(&boundary::is_in_conductor,
									      bound,
									      std::placeholders::_1,
									      std::placeholders::_2,
									      std::placeholders::_3));
  is_dirichlet = new field <bool> (x_axis, y_axis, z_axis, is_dirichlet_func);

  // is this point inside of a von neumann BC?
  std::function <bool (double, double, double)> is_vn_func;
  is_vn_func = std::function<bool (double, double, double)> (std::bind(&boundary::is_in_VN,
								       bound,
								       std::placeholders::_1,
								       std::placeholders::_2,
								       std::placeholders::_3));
  is_von_neumann = new field <bool> (x_axis, y_axis, z_axis, is_vn_func);

  is_boundary = new field <bool> (x_axis, y_axis, z_axis, false);
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	is_boundary -> set(i, j, k,
			   ( is_dirichlet -> get(i, j, k) )
			   or ( is_von_neumann -> get(i, j, k) ));
      }
    }
  }
  
  std::function <double (double, double, double)> vn_E_func;
  vn_E_func = std::function<double (double, double, double)> (std::bind(&boundary::Efield,
									bound,
									std::placeholders::_1,
									std::placeholders::_2,
									std::placeholders::_3));
  von_neumann_dV = new field <double> (x_axis, y_axis, z_axis, vn_E_func);
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	von_neumann_dV -> set(i, j, k, (von_neumann_dV -> get(i, j, k))*spacing);
      }
    }
  }
}

void solver::initialize_permittivity_and_conductivity()
{
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

void solver::solve_static()
{
  std::cout << "#\n"
	    << "# Iteration:"
	    << '\t'
	    << "Stepwise Difference Norm:"
	    << std::endl;
  for ( int iter = 0; iter < nIter; iter++ ) {
    relax(1);

    set_VN();

    if ( sig_int ) {
      // check for ^C
      std::cout << std::endl
		<< "# Recieved SIGTERM! Saving and quitting..."
		<< std::endl;
      break;
    }

    int done = report(iter);
    if ( done ) break;
    
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
  
  // field <double> * VNtemp = new field <double> (x_axis, y_axis, y_axis, constant(0.));
  fill_empty_VN();

  std::cout << "# final stepwise difference: "
	    << squared_sum(resid, is_boundary)
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

void solver::set_VN()
{
  // set the value of the potential where it is von Neumann
  for ( int i = 0; i < nPointsX; i++ ) {
    for ( int j = 0; j < nPointsY; j++ ) {
      for ( int k = 0; k < nPointsZ; k++ ) {
	if ( not is_boundary -> get(i, j, k) ) {
	  if ( is_von_neumann -> get(i+1, j, k) ) {
	    potential -> set(i+1, j, k,
			     potential -> get(i, j, k) - von_neumann_dV -> get(i+1, j, k));
	  }
	  if ( is_von_neumann -> get(i-1, j, k) ) {
	    potential -> set(i-1, j, k,
			     potential -> get(i, j, k) + von_neumann_dV -> get(i-1, j, k));
	  }
	  if ( is_von_neumann -> get(i, j+1, k) ) {
	    potential -> set(i, j+1, k,
			     potential -> get(i, j, k) - von_neumann_dV -> get(i, j+1, k));
	  }
	  if ( is_von_neumann -> get(i, j-1, k) ) {
	    potential -> set(i, j-1, k,
			     potential -> get(i, j, k) + von_neumann_dV -> get(i, j-1, k));
	  }
	  if ( is_von_neumann -> get(i, j, k+1) ) {
	    potential -> set(i, j, k+1,
			     potential -> get(i, j, k) - von_neumann_dV -> get(i, j, k+1));
	  }
	  if ( is_von_neumann -> get(i, j, k-1) ) {
	    potential -> set(i, j, k-1,
			     potential -> get(i, j, k) + von_neumann_dV -> get(i, j, k-1));
	  }
	}
      }
    }
  }
}

void solver::fill_empty_VN()
{
  // need to do this twice
  // for the voxels in corners
  for ( int iter = 0; iter < 2; iter++ ) {
    for ( int i = 0; i < potential -> xSize; i++ ) {
      for ( int j = 0; j < potential -> ySize; j++ ) {
	for ( int k = 0; k < potential -> zSize; k++ ) {
	  if ( potential -> get(i, j, k) == 0xdeadbeef ) {
	    // this is a piece of von neumann boundary that never
	    // got set, since it doesn't neighbor any non-boundary
	    // voxels
	    double neighborSum = 0;
	    int nNeighbors = 0;
	    if ( i != (potential -> xSize - 1) ) {
	      if ( ( is_von_neumann -> get(i+1, j, k) )
		   and ( potential -> get(i+1, j, k) != 0xdeadbeef ) ) {
		neighborSum += potential -> get(i+1, j, k);
		nNeighbors++;
	      }
	    }
	    if ( i != 0 ) {
	      if ( ( is_von_neumann -> get(i-1, j, k) )
		   and ( potential -> get(i-1, j, k) != 0xdeadbeef ) ) {
		neighborSum += potential -> get(i-1, j, k);
		nNeighbors++;
	      }
	    }
	    if ( j != (potential -> ySize - 1) ) {
	      if ( ( is_von_neumann -> get(i, j+1, k) )
		   and ( potential -> get(i, j+1, k) != 0xdeadbeef )) {
		neighborSum += potential -> get(i, j+1, k);
		nNeighbors++;
	      }
	    }
	    if ( j != 0 ) {
	      if ( ( is_von_neumann -> get(i, j-1, k) )
		   and ( potential -> get(i, j-1, k) != 0xdeadbeef )) {
		neighborSum += potential -> get(i, j-1, k);
		nNeighbors++;
	      }
	    }
	    if ( k != (potential -> zSize - 1) ) {
	      if ( ( is_von_neumann -> get(i, j, k+1) )
		   and ( potential -> get(i, j, k+1) != 0xdeadbeef )) {
		neighborSum += potential -> get(i, j, k+1);
		nNeighbors++;
	      }
	    }
	    if ( k != 0 ) {
	      if ( ( is_von_neumann -> get(i, j, k-1) )
		   and ( potential -> get(i, j, k-1) != 0xdeadbeef )) {
		neighborSum += potential -> get(i, j, k-1);
		nNeighbors++;
	      }
	    }
	    if ( nNeighbors ) {
	      von_neumann_dV -> set(i, j, k, neighborSum/nNeighbors);
	    }
	    else {
	      von_neumann_dV -> set(i, j, k, 0xdeadbeef);
	    }
	  }
	}
      }
    }
    for ( int i = 0; i < potential -> xSize; i++ ) {
      for ( int j = 0; j < potential -> ySize; j++ ) {
	for ( int k = 0; k < potential -> zSize; k++ ) {
	  if ( potential -> get(i, j, k) == 0xdeadbeef ) {
	    potential -> set(i, j, k, von_neumann_dV -> get(i, j, k));
	  }
	}
      }
    }
  }
}

void solver::relax(int nThreads = 2)
{
  omp_set_dynamic(0);
  omp_set_num_threads(4);

  #pragma omp parallel
  {

    // std::cout << "# threads: " << omp_get_num_threads() << std::endl;
    #pragma omp for
    for ( int i = 0; i < potential -> xSize; i++ ) {
      // std::cout << omp_get_num_threads() << '\t'
      // 		<< omp_get_thread_num() << '\t'
      // 		<< i << std::endl;
      for ( int j = 0; j < potential -> ySize; j++ ) {
	for ( int k = 0; k < potential -> zSize; k++ ) {
	  if ( not  is_boundary -> get(i, j, k) ) {
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
  }
}

int solver::report(int iter)
{
  if ( iter % iterFreq == 0 ) {
    // print out the squared_diff between current iteration and previous iteration
    double sqsm = squared_sum(resid, is_boundary);
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
      return 1;
    }
  }
  return 0;
}

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);

  // omp_set_num_threads(4);
  
  // sig handling
  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = term;
  sigaction(SIGINT, &action, NULL);
  
  boundary * detector = new boundary (geom);
  solver thisSolver (detector, N, initSpacing);
  thisSolver.sigVal -> print_to_file("cond_map.dat");
  thisSolver.is_dirichlet -> print_to_file("dirichlet_map.dat");
  thisSolver.is_von_neumann -> print_to_file("VN_map.dat");
  thisSolver.epVal -> print_to_file("perm_map.dat");
  thisSolver.potential -> print_to_file("initial.dat");
  
  // thisSolver.solve_charge();
  thisSolver.solve_static();
  // for ( int order = 0; order < 3; order++ ) {
  //   thisSolver.upscale(2);
  //   thisSolver.solve_static();
  // }
  thisSolver.potential -> print_to_file(outFileName);
  thisSolver.Q -> print_to_file("Q.dat");
    
  return 0;
}
