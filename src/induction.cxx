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
std::string geom = "geometries/bulkPix.json";

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
  boundary detector (geom);

  const int nPads = 4;
  // How finely to sample the drift paths
  // 1 is just the central (head-on) path
  const int fillFac = 3;
  const int pathsPerPad = (2*(fillFac - 1) + 1);
  const double padPitch = 0.4;

  // Set up a grid so that we can reuse drift paths by translation and reflection
  // We only need to simulate drift above one quadrant of one pad
  // This grid must be within the potential field
  const int nPaths = fillFac;
  const double endPoint = (fillFac - 1)/float(2*fillFac - 1)*padPitch;

  // The larger grid where induction is calculated
  // This grid must be within the weighting field
  const int nPathsInduct = 2*nPads*fillFac - nPads - fillFac + 1;
  const double endPointInduct = (nPads - 1 + (fillFac - 1)/float(2*fillFac - 1))*padPitch;

  // std::thread * threads [nPaths*nPaths];
  path * trajectories [nPathsInduct*nPathsInduct];
  bool trajIsSet [nPathsInduct*nPathsInduct] = {false};
  
  std::vector <double> xi_space = linspace(0, endPoint, nPaths);
  std::vector <double> yi_space = linspace(0, endPoint, nPaths);
  double zi = 0.4;
  
  for ( int i = 0; i < xi_space.size(); i++ ) {
    double xi = xi_space[i];
    for ( int j = 0; j < yi_space.size(); j++ ) {
      double yi = yi_space[j];
      drift_path(std::vector <double> {xi, yi, zi},
		 potential,
		 detector,
		 std::ref(trajectories[i*nPathsInduct + j]));
      trajIsSet[i*nPathsInduct + j] = true;
      // threads[i*nPaths + j] = new std::thread ( drift_path,
      // 						std::vector <double> {xi, yi, zi},
      // 						potential,
      // 						detector,
      // 						std::ref(trajectories[i*nPaths + j]) );
    }
  }

  std::vector <double> xiInduct_space = linspace(0, endPointInduct, nPathsInduct);
  std::vector <double> yiInduct_space = linspace(0, endPointInduct, nPathsInduct);

  for ( int i = 0; i < xiInduct_space.size(); i++ ) {
    double xi = xiInduct_space[i];
    for ( int j = 0; j < yiInduct_space.size(); j++ ) {
      double yi = yiInduct_space[j];

      // Find out which pad this starting point corresponds to 
      if ( not trajIsSet[i*nPathsInduct + j] ) {
	int iOffset = 0;
	while ( i - iOffset*pathsPerPad >= fillFac ) {
	  iOffset++;
	}

	int jOffset = 0;
	while ( j - jOffset*pathsPerPad >= fillFac ) {
	  jOffset++;
	}

	// Is this point on the left (lower) or right (upper) side of the pad?
	// If left (lower), we need to reflect the precalculated path
	bool iReflect = (i - iOffset*pathsPerPad < 0);
	bool jReflect = (j - jOffset*pathsPerPad < 0);

	// The index of the corresponding precalculated path
	int iPreCalc = abs(i - iOffset*pathsPerPad);
	int jPreCalc = abs(j - jOffset*pathsPerPad);

	path * thisTraj = trajectories[iPreCalc*nPathsInduct + jPreCalc] -> copy();
	if ( iReflect ) thisTraj -> reflectX();
	if ( jReflect ) thisTraj -> reflectY();
	thisTraj -> shift({iOffset*padPitch, jOffset*padPitch, 0});

	trajectories[i*nPathsInduct + j] = thisTraj;
	
	// std::cout << i << ", " << xi << ", " << iOffset << ", " << iReflect << ", " << iPreCalc << '\t'
	// 	  << j << ", " << yi << '\t' << jOffset << ", " << jReflect << ", " << jPreCalc
	// 	  << std::endl;
      }
    }
    std::ostringstream stringStream1;
    stringStream1 << "drift_from_" << xi << ".dat";
    trajectories[i*nPathsInduct] -> print_to_file(stringStream1.str());

    std::ostringstream stringStream2;
    stringStream2 << "wf_from_" << xi << ".dat";
    std::vector <double> pSeries = ramo_induction(e,
						  trajectories[i*nPathsInduct],
						  weight);
    std::ofstream outFile (stringStream2.str().c_str());
    for ( int k = 0; k < pSeries.size(); k++ ) {
      outFile << pSeries[k] << '\n';
    }
    outFile.close();
  }
  
  // for ( int i = 0; i < xi_space.size(); i++ ) {
  //   double xi = xi_space[i];
  //   for ( int j = 0; j <= i; j++ ) {
  //     double yi = yi_space[j];
  //     threads[i*nPaths + j] -> join();
  //     // std::vector <double> pSeries = ramo_induction(e,
  //     // 						    trajectories[i*nPaths + j],
  //     // 						    weight);
  //     std::cout << "done the thing here " << i << '\t' << j << std::endl;
  //   }
  // }
  
  // int Nsegments = 50;
  // double t0 = 0;
  // double dx = 0.01;
  // path * driftPaths [Nsegments];
  // for ( int i = 0; i < Nsegments; i++ ) {
  //   pad * centerPad = new pad(150, std::vector <double> {0,0,0});

  //   double xi = i*dx;
  //   std::vector <double> startingPos = {xi, 0, 0};
  //   std::cout << startingPos[0] << std::endl;
    
  //   driftPaths[i] = new path (dt);
  //   drift_path(startingPos, potential, bound, driftPaths[i]);

  //   centerPad -> add_response(e, t0, startingPos, driftPaths[i],
  // 			      potential, weight, bound);

  //   std::stringstream responseFile;
  //   responseFile << "responses/" << xi << ".dat";
  //   centerPad -> print_current_to_file(responseFile.str());

  //   std::stringstream driftPathFile;
  //   driftPathFile << "paths/" << xi << ".dat";
  //   driftPaths[i] -> print_to_file(driftPathFile.str());
    
  //   delete centerPad;
  // }

  return 0;
}
