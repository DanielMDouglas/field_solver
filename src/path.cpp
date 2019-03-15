#include <fstream>
#include <iostream>

#include "path.h"

path::path(double deltaT)
{
  dt = deltaT;
}

path::path(std::string filename)
{
  std::ifstream inFile (filename.c_str());

  std::string entry;
  double value = 0;
  int nSteps = 0;

  getline(inFile, entry, ',');
  dt = std::stod(entry);
  
  getline(inFile, entry, ',');
  arrivalTime = std::stod(entry);
  
  getline(inFile, entry, '\n');
  fate = entry;

  while ( getline(inFile, entry, ',') ) {
  
    std::vector <double> step;

    getline(inFile, entry, ',');
    step.push_back(std::stod(entry));

    getline(inFile, entry, ',');
    step.push_back(std::stod(entry));

    getline(inFile, entry, '\n');
    step.push_back(std::stod(entry));

    steps.push_back(step);
  }
}

void path::print_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  outFile << dt << ','
	  << arrivalTime << ','
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
