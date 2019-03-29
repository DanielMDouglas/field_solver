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

std::string arrTimeFileName = "none";
std::string pSeriesFileName = "none";

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
    if ( optValue.str() == "-oa" ) {
      argValue >> arrTimeFileName;
    }
    if ( optValue.str() == "-op" ) {
      argValue >> pSeriesFileName;
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
	    << "arrival time output: " << arrTimeFileName << '\n'
	    << "pulse series output: " << pSeriesFileName << '\n'
	    << std::endl;
}

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);
  
  scalarField <double> * potential = new scalarField <double> (potFileName);
  scalarField <double> * weight = new scalarField <double> (weightFileName);
  boundary bound (geom);

  pad * padList [85];
  int nPads = 0;
  
  for ( int i = 0; i < 17; i++ ) {
    double x = -3.2 + 0.4*i;
    for ( int j = 0; j < 5; j++ ) {
      double y = -0.8 + 0.4*j;
      std::cout << x << '\t' << y << std::endl;
      padList[nPads] = new pad(150, std::vector <double> {x, y, 0});
      nPads++;
    }
  }
  
  // pad * centerPad = new pad(150, std::vector <double> {0,0,0});
  // pad * adjPad = new pad(150, std::vector <double> {0.4, 0, 0});
  // pad * diagPad = new pad(150, std::vector <double> {0.4, 0.4, 0});

  std::ifstream inFile ("track_sample_small.dat");
  std::string entry;
  double value = 0;
  
  while ( getline(inFile, entry, '\t') ) {
    std::vector <double> init_pos;
    init_pos.push_back(std::stod(entry));

    getline(inFile, entry, '\t');
    init_pos.push_back(std::stod(entry));

    getline(inFile, entry, '\t');
    init_pos.push_back(std::stod(entry));

    getline(inFile, entry, '\n');
    double t0 = std::stod(entry);

    std::cout << "New drift from "
	      << "("  << init_pos[0]
	      << ", " << init_pos[1]
	      << ", " << init_pos[2]
	      << ")"  << std::endl;

    std::vector <double> driftPos = {fmod(init_pos[0], 0.4),
				     fmod(init_pos[1], 0.4),
				     init_pos[2]};

    std::cout << "Calculating path from "
	      << "("  << driftPos[0]
	      << ", " << driftPos[1]
	      << ", " << driftPos[2]
	      << ")"  << std::endl;
    
    path * driftPath;
    drift_path(driftPos, potential, bound, driftPath);

    driftPath -> shift(init_pos - driftPos);
    
    for ( pad * thisPad: padList ) {
      thisPad -> add_response(-e, t0, init_pos, driftPath,
			      potential, weight, bound);
    }

    delete driftPath;
  }

  inFile.close();

  for ( pad * thisPad: padList) {
    std::ostringstream stringStream;
    stringStream << "responses_small/"
		 << thisPad -> center[0] << "_"
		 << thisPad -> center[1] << "_response.dat";
    thisPad -> print_to_file(stringStream.str());
  }
    
  // for ( int i = 0; i < 100000; i++ ) {
  //   double * t;
  //   std::vector <double> * samplePos;
  //   sample(std::vector <double> {0, 0, 10.}, t, samplePos);

  //   // std::cout << &t << '\t' << &samplePos[0] << std::endl;
  // }
  
  return 0;
}
