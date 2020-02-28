#include <fstream>
#include <iostream>

#include "path.h"

path::path(double deltaT)
{
  dt = deltaT;
  // steps = std::vector <std::vector <double>> (1e6, std::vector <double> (3));
}

path::path(std::string filename)
{
  std::ifstream inFile (filename.c_str());

  std::string entry;
  double value = 0;
  int nSteps = 0;

  // first thing is dt
  getline(inFile, entry, ',');
  dt = std::stod(entry);

  // final time
  getline(inFile, entry, ',');
  arrivalTime = std::stod(entry);

  // fate
  getline(inFile, entry, '\n');
  fate = entry;

  // each line begins with time
  while ( getline(inFile, entry, ',') ) {
  
    std::vector <double> step;

    // x position
    getline(inFile, entry, ',');
    step.push_back(std::stod(entry));

    // y position
    getline(inFile, entry, ',');
    step.push_back(std::stod(entry));

    // z position
    getline(inFile, entry, ',');
    step.push_back(std::stod(entry));

    pos.push_back(step);

    std::vector <double> v;
	
    // x velocity
    getline(inFile, entry, ',');
    v.push_back(std::stod(entry));

    // y velocity
    getline(inFile, entry, ',');
    v.push_back(std::stod(entry));

    // z velocity
    getline(inFile, entry, '\n');
    v.push_back(std::stod(entry));

    vel.push_back(v);
  }
}

path * path::copy()
{
  path * newPath = new path(dt);
  newPath -> pos = pos;
  newPath -> vel = vel;
  newPath -> E = E;
  newPath -> fate = fate;
  newPath -> arrivalTime = arrivalTime;
  
  return newPath;
}

void path::reflectX()
{
  for ( uint i = 0; i < pos.size(); i++ ) {
    pos[i][0] = -pos[i][0];
    vel[i][0] = -vel[i][0];
    E[i][0] = -E[i][0];
  }

}

void path::reflectY()
{
  for ( uint i = 0; i < pos.size(); i++ ) {
    pos[i][1] = -pos[i][1];
    vel[i][1] = -vel[i][1];
    E[i][1] = -E[i][1];
  }
}

void path::reflectZ()
{
  for ( uint i = 0; i < pos.size(); i++ ) {
    pos[i][2] = -pos[i][2];
    vel[i][2] = -vel[i][2];
    E[i][2] = -E[i][2];
  }
}

void path::shift(std::vector <double> offset)
{
  for ( uint i = 0; i < pos.size(); i++ ) {
    pos[i] = pos[i] + offset;
  }
}

void path::print_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  outFile << dt << ','
	  << arrivalTime << ','
	  << fate << '\n';
  for ( uint i = 0; i < pos.size(); i++ ) {
    double t = i*dt;
    outFile << t << ','
	    << pos[i][0] << ','
	    << pos[i][1] << ','
	    << pos[i][2] << ','
	    << vel[i][0] << ','
	    << vel[i][1] << ','
	    << vel[i][2] << ','      
	    << E[i][0] << ','
	    << E[i][1] << ','
	    << E[i][2] << '\n';      
  }
  outFile.close();
}
