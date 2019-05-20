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
  
  // pad * centerPad = new pad(150, std::vector <double> {0,0,0});
  // pad * adjPad = new pad(150, std::vector <double> {0.4, 0, 0});
  // pad * diagPad = new pad(150, std::vector <double> {0.4, 0.4, 0});

  int nPads = 16*16;
  pad * padList [nPads];
  int i = 0;
  double pitch = 0.4;
  for ( double x = 7.6; x >= 7.6 - 16*pitch; x -= pitch ) {
    for ( double y = 5.2; y >= 5.2 - 16*pitch; y -= pitch ) {
      padList[i] = new pad(10, std::vector <double> {x, y, 0});
      i++;
    }
  }
      
  std::vector <double> trackOrigPos = {5, 4, 1.2};
  std::vector <double> trackDir = {-0.707, -0.707, 0};
  double trackLen = 6;
  int Nsegments = 20;
  double t0 = 0;
  path * driftPaths [Nsegments];
  for ( int i = 0; i < Nsegments; i++ ) {
    double distance = trackLen*(float(i)/Nsegments);
    std::vector <double> startingPos = trackOrigPos + trackDir*distance;
    std::cout << startingPos[0] << std::endl;
    std::vector <double> relativePos = {fmod(startingPos[0] + 0.2, 0.4) - 0.2,
					fmod(startingPos[1] + 0.2, 0.4) - 0.2,
					startingPos[2]};
    
    driftPaths[i] = new path (dt);
    drift_path(relativePos , potential, bound, driftPaths[i]);
    driftPaths[i] -> shift(startingPos - relativePos);
    for ( pad * thisPad: padList ) {
      thisPad -> add_response(-20*e, t0, startingPos, driftPaths[i],
			      potential, weight, bound);
    }
  }

  padList[141] -> print_current_to_file("currentSeries.dat");
  
  for ( pad * thisPad: padList ) {
    thisPad -> calculate_output();
    thisPad -> print_output_to_file(outFileName);
  }

  return 0;
}
