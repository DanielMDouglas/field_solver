#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "vector_field.h"

vectorField::vectorField(scalarField phi, std::string type = "grad")
{
  xSize = phi.xSize;
  ySize = phi.ySize;
  zSize = phi.zSize;

  x_space = phi.x_space;
  y_space = phi.y_space;
  z_space = phi.z_space;

  values = std::vector<std::vector <double>> (xSize*ySize*zSize, std::vector <double> (3, 0));

  // assume constant grid spacing
  double dx = x_space[1] - x_space[0];
  double dy = y_space[1] - y_space[0];
  double dz = z_space[1] - z_space[0];
  
  if ( type == "grad" ) {
    for ( int i = 0; i < xSize; i++ ) {
      for ( int j = 0; j < ySize; j++ ) {
  	for ( int k = 0; k < zSize; k++ ) {
	  double ddx = (phi.get(i+1, j, k) - phi.get(i-1, j, k))/(2*dx);
	  double ddy = (phi.get(i, j+1, k) - phi.get(i, j-1, k))/(2*dy);
	  double ddz = (phi.get(i, j, k+1) - phi.get(i, j, k-1))/(2*dz);
	  std::vector <double> diff = {ddx, ddy, ddz};
	  set(i, j, k, diff);
  	}
      }
    }
  }
}

void vectorField::set(int i, int j, int k, std::vector <double> value)
{
  values[ySize*zSize*i + zSize*j + k] = value;
}

std::vector <double> vectorField::get(int i, int j, int k)
{
  return values[ySize*zSize*i + zSize*j + k];
}

void vectorField::print_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      for ( int k = 0; k < zSize; k++ ) {
	std::vector <double> thisValue = get(i, j, k);
	outFile << x_space[i] << ','
		<< y_space[j] << ','
		<< z_space[k] << ','
		<< thisValue[0] << ','
		<< thisValue[1] << ','
		<< thisValue[2] << '\n';
      }
    }
  }
  outFile.close();
}
