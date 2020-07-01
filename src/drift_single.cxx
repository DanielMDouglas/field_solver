#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <math.h>
#include <sstream>

#include "drift_single.h"

std::string fieldFileName = "none";
std::string geom = "geometries/bulkPix.json";
std::string outFileName = "drift.dat";
double xi = 0;
double yi = 0;
double zi = 0.45;

void handleOpts(int argc, const char ** argv)
{
  argParser parser;

  parser.add_option("-f", [](arg_t ss) {*ss >> fieldFileName;});
  parser.add_option("-g", [](arg_t ss) {*ss >> geom;});
  parser.add_option("-o", [](arg_t ss) {*ss >> outFileName;});
  parser.add_option("-x", [](arg_t ss) {*ss >> xi;});
  parser.add_option("-y", [](arg_t ss) {*ss >> yi;});
  parser.add_option("-z", [](arg_t ss) {*ss >> zi;});

  parser.add_flag("-h", sayUsage);

  parser.parse(argc, argv);
  
  if ( fieldFileName == "none" ) {
    std::cout << "Need a field file!" << std::endl;
    exit(1);
  }
}

void saySettings()
{
  std::cout << "#####################################" << '\n'
	    << "# Using arguments: \n"
	    << "# xi:               " << xi << '\n'
	    << "# yi:               " << yi << '\n'
	    << "# zi:               " << zi << '\n'
	    << "# outFile:          " << outFileName << '\n'
	    << "# potential field:  " << fieldFileName << '\n'
	    << "# geom:             " << geom << '\n'
	    << "#####################################" << std::endl;
}

void sayUsage(arg_t ss)
{
  std::cout << "Usage: ./drift_single [OPTIONS]" << std::endl
	    << "-x" << '\t' << "Initial x position, default 0" << std::endl
	    << "-y" << '\t' << "Initial y position, default 0" << std::endl
	    << "-z" << '\t' << "Initial z position, default 0.45" << std::endl
	    << "-o" << '\t' << "output, default: drift.dat" << std::endl
	    << "-f" << '\t' << "electric potential field, default: none" << std::endl
	    << "-g" << '\t' << "detector geometry, default: geometries/bulkPix.json" << std::endl;
}  

int main(int argc, const char ** argv)
{
  handleOpts(argc, argv);
  
  field <double> * potential = new field <double> (fieldFileName);
  boundary detector (geom);

  path* thisTraj;
  drift_path(std::vector <double> {xi, yi, zi},
	     potential,
	     detector,
	     thisTraj);

  std::ostringstream stringStream;
  stringStream << "drift_from_" << xi << ".dat";
  thisTraj -> print_to_file(stringStream.str());
  
  return 0;
}
