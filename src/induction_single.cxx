#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <math.h>

#include "induction.h"
#include "physics.h"
#include "pad.h"

std::string weightFileName = "none";
std::string driftFileName = "none";

std::string outFileName = "current.dat";

void handleOpts(int argc, const char ** argv)
{
  argParser parser;

  parser.add_option("-d", [](arg_t ss) {*ss >> driftFileName;});
  parser.add_option("-w", [](arg_t ss) {*ss >> weightFileName;});
  parser.add_option("-o", [](arg_t ss) {*ss >> outFileName;});

  parser.add_flag("-h", [](arg_t ss) {sayUsage(ss);});
    
  parser.parse(argc, argv);

  if ( driftFileName == "none" ) {
    std::cout << "Need a drift path file!" << std::endl;
    exit(1);
  }
  if ( weightFileName == "none" ) {
    std::cout << "Need a weighting field file!" << std::endl;
    exit(1);
  }

  saySettings();
}

void saySettings()
{
  std::cout << "#####################################" << '\n'
	    << "# Using arguments: \n"
	    << "# drift path:          " << driftFileName << '\n'
	    << "# weighting field:     " << weightFileName << '\n'
	    << "# output:              " << outFileName << '\n'
	    << "#####################################" << '\n'
	    << std::endl;
}

void sayUsage(arg_t ss)
{
  std::cout << "Usage: ./induction_single [OPTIONS]" << std::endl
	    << "-d" << '\t' << "potential field, default: none" << std::endl
	    << "-w" << '\t' << "weighting field, default: none" << std::endl
	    << "-o" << '\t' << "output, default: current.dat" << std::endl
	    << "-h" << '\t' << "display this help and exit" << std::endl;
}

int main(int argc, const char ** argv)
{
  handleOpts(argc, argv);
  
  field <double> * weight = new field <double> (weightFileName);
  path * driftPath = new path(driftFileName);
  
  std::vector <double> pSeries = ramo_induction(-e,
						driftPath,
						weight);
  std::ofstream outFile (outFileName);
  for ( int k = 0; k < pSeries.size(); k++ ) {
    outFile << k*dt << '\t' << pSeries[k] << '\n';
  }
  outFile.close();

  return 0;
}
