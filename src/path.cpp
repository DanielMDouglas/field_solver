#include <fstream>

#include "path.h"

path::path(double deltaT)
{
  dt = deltaT;
}

void path::print_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  outFile << arrivalTime << ','
	  << fate << '\n';
  for ( uint i = 0; i < steps.size(); i++ ) {
    double t = i*dt;
    outFile << t << ','
	    << steps[i][0] << ','
	    << steps[i][1] << ','
	    << steps[i][2] << '\n';      
  }
  outFile.close();
}
