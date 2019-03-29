#include <iostream>
#include <string>
#include <sstream>

#include "efield.h"

std::string inFileName = "none";
std::string outFileBase = "none";

void handleOpts(int argc, char const * argv[])
{
  int opt = 0;
  while ( opt < argc ) {
    std::stringstream optValue;
    std::stringstream argValue;
    optValue << argv[++opt];
    argValue << argv[++opt];
    
    if ( optValue.str() == "-i" ) {
      argValue >> inFileName;
    }
    if ( optValue.str() == "-o" ) {
      argValue >> outFileBase;
    }
  }

  if ( inFileName == "none" ) {
    std::cout << "Need an input file!" << std::endl;
    exit(1);
  }
  if ( outFileBase == "none" ) {
    std::cout << "Need an output file!" << std::endl;
    exit(1);
  }
  
  std::cout << "Using arguments: \n"
	    << "inFile:  " << inFileName << '\n'
	    << "outFileBase: " << outFileBase << std::endl;
}

int main(int argc, char const * argv[])
{
  handleOpts(argc, argv);

  scalarField <double> potential (inFileName);

  scalarField <double> Ex (potential.x_space,
			   potential.y_space,
			   potential.z_space,
			   0);
  scalarField <double> Ey (potential.x_space,
			   potential.y_space,
			   potential.z_space,
			   0);
  scalarField <double> Ez (potential.x_space,
			   potential.y_space,
			   potential.z_space,
			   0);
  scalarField <double> Emag (potential.x_space,
			     potential.y_space,
			     potential.z_space,
			     0);

  for ( int i = 0; i < potential.xSize; i++ ) {
    for ( int j = 0; j < potential.ySize; j++ ) {
      for ( int k = 0; k < potential.zSize; k++ ) {
	std::vector <double> pos = {potential.x_space[i],
				    potential.y_space[j],
				    potential.z_space[k]};
	std::vector <double> Efield = E(pos, &potential);
	Ex.set(i, j, k, Efield[0]);
	Ey.set(i, j, k, Efield[1]);
	Ez.set(i, j, k, Efield[2]);
	Emag.set(i, j, k, mag(Efield));
      }
    }
  }

  std::stringstream outFile;
  outFile << outFileBase << "_Ex.dat";
  Ex.print_to_file(outFile.str());

  outFile.str("");
  outFile << outFileBase << "_Ey.dat";
  Ey.print_to_file(outFile.str());

  outFile.str("");
  outFile << outFileBase << "_Ez.dat";
  Ez.print_to_file(outFile.str());

  outFile.str("");
  outFile << outFileBase << "_Emag.dat";
  Emag.print_to_file(outFile.str());
  
  return 0;
}
