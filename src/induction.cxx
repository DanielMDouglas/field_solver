#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <math.h>

#include "induction.h"
#include "physics.h"
#include "pad.h"

std::string potFileName = "none";
std::string weightFileName = "none";
std::string geom = "bulkPix";

std::string outFileName = "event.dat";

void handleOpts(int argc, char const * argv[])
{
  int opt = 0;
  while ( opt < argc ) {
    std::stringstream optValue;
    std::stringstream argValue;
    optValue << argv[++opt];
    argValue << argv[++opt];
    
    if ( optValue.str() == "-p" ) {
      argValue >> potFileName;
    }
    if ( optValue.str() == "-w" ) {
      argValue >> weightFileName;
    }
    if ( optValue.str() == "-g" ) {
      argValue >> geom;
    }
    if ( optValue.str() == "-o" ) {
      argValue >> outFileName;
    }
  }

  if ( potFileName == "none" ) {
    std::cout << "Need a potential field file!" << std::endl;
    exit(1);
  }
  if ( weightFileName == "none" ) {
    std::cout << "Need a weighting field file!" << std::endl;
    exit(1);
  }

  std::cout << "Using arguments: \n"
	    << "potential field:     " << potFileName << '\n'
	    << "weighting field:     " << weightFileName << '\n'
	    << "geometry:            " << geom << '\n'
	    << "output:              " << outFileName << '\n'
	    << std::endl;
}

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);
  
  field <double> * potential = new field <double> (potFileName);
  field <double> * weight = new field <double> (weightFileName);
  boundary bound (geom);
  
  int Nsegments = 50;
  double t0 = 0;
  double dx = 0.01;
  path * driftPaths [Nsegments];
  for ( int i = 0; i < Nsegments; i++ ) {
    pad * centerPad = new pad(150, std::vector <double> {0,0,0});

    double xi = i*dx;
    std::vector <double> startingPos = {xi, 0, 0};
    std::cout << startingPos[0] << std::endl;
    
    driftPaths[i] = new path (dt);
    drift_path(startingPos, potential, bound, driftPaths[i]);

    centerPad -> add_response(e, t0, startingPos, driftPaths[i],
			      potential, weight, bound);

    std::stringstream responseFile;
    responseFile << "responses/" << xi << ".dat";
    centerPad -> print_current_to_file(responseFile.str());

    std::stringstream driftPathFile;
    driftPathFile << "paths/" << xi << ".dat";
    driftPaths[i] -> print_to_file(driftPathFile.str());
    
    delete centerPad;
  }

  return 0;
}
