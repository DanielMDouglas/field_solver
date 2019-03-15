#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "induction.h"
#include "physics.h"

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
  
  scalarField * potential = new scalarField(potFileName);
  scalarField * weight = new scalarField(weightFileName);
  boundary bound (geom);

  std::ofstream arrTimeOutFile (arrTimeFileName);
  std::ofstream pSeriesOutFile (pSeriesFileName);
  
  double ds = 0.01;
  for ( double xi = 0; xi <= 0.2; xi += ds ) {
    for ( double yi = 0; yi <= xi; yi += ds ) {
      std::vector <double> init_pos = {xi, yi, 2.};
      path driftPath = drift_path(init_pos, potential, bound);
      if ( xi <= 0.2 ) {
	arrTimeOutFile << xi << ',' << yi << ',' << driftPath.steps.size()*driftPath.dt << '\n';
      }
      std::vector <double> pSeries = ramo_induction(driftPath, weight);
      double maxAbs = 0;
      for ( double val: pSeries ) {
	if ( abs(val) > maxAbs ) {
	  maxAbs = abs(val);
	}
      }
      pSeriesOutFile << xi << ',' << yi << ',' << maxAbs << '\n';
    }
  }

  arrTimeOutFile.close();
  pSeriesOutFile.close();
  
  return 0;
}
