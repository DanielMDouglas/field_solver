#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <math.h>

#include "scalar_field.h"

scalarField::scalarField(std::vector <double> x,
			 std::vector <double> y,
			 std::vector <double> z,
			 double init_val = 0)
{
  xSize = x.size();
  ySize = y.size();
  zSize = z.size();
  
  x_space = x;
  y_space = y;
  z_space = z;

  values = std::vector< double > (xSize*ySize*zSize, init_val);
}

scalarField::scalarField(std::vector <double> x,
			 std::vector <double> y,
			 std::vector <double> z,
			 std::function <double (double, double, double)> f)
{
  xSize = x.size();
  ySize = y.size();
  zSize = z.size();

  x_space = x;
  y_space = y;
  z_space = z;

  values = std::vector <double> (xSize*ySize*zSize);

  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      for ( int k = 0; k < zSize; k++ ) {
	set(i, j, k, f(x[i], y[j], z[k]));
      }
    }
  }
}

scalarField::scalarField(boundary bound, int Nx, int Ny, int Nz, std::string type = "val")
{
  xSize = Nx;
  ySize = Ny;
  zSize = Nz;

  x_space = linspace(bound.Xmin, bound.Xmax, Nx);
  y_space = linspace(bound.Ymin, bound.Ymax, Ny);
  z_space = linspace(bound.Zmin, bound.Zmax, Nz);

  values = std::vector< double > (xSize*ySize*zSize);

  if ( type == "val" ) {
    for ( int i = 0; i < xSize; i++ ) {
      for ( int j = 0; j < ySize; j++ ) {
	for ( int k = 0; k < zSize; k++ ) {
	  set(i, j, k,
	      bound.boundary_value(x_space[i],
				   y_space[j],
				   z_space[k]));
	}
      }
    }
  }
  else if ( type == "bool" ) {
    for ( int i = 0; i < xSize; i++ ) {
      for ( int j = 0; j < ySize; j++ ) {
	for ( int k = 0; k < zSize; k++ ) {
	  if ( bound.is_in_boundary(x_space[i],
				    y_space[j],
				    z_space[k]) ) {
	    set(i, j, k, 1);
	  }
	  else {
	    set(i, j, k, 0);
	  }
	}
      }
    }
  }  
}

scalarField::scalarField(std::string filename)
{
  std::ifstream inFile (filename.c_str());

  std::string entry;
  double value;
  int nLines = 0;

  while (getline(inFile, entry, ',')) {
    value = std::stod(entry);
    if ( std::find ( x_space.begin(), x_space.end(), value ) == x_space.end() ) {
      x_space.push_back(value);
    }
    
    getline(inFile, entry, ',');
    value = std::stod(entry);
    if ( std::find ( y_space.begin(), y_space.end(), value ) == y_space.end() ) {
      y_space.push_back(value);
    }

    getline(inFile, entry, ',');
    value = std::stod(entry);
    if ( std::find ( z_space.begin(), z_space.end(), value ) == z_space.end() ) {
      z_space.push_back(value);
    }

    getline(inFile, entry, '\n');
    values.push_back(std::stod(entry));
    
    nLines++;
  }

  inFile.close();
  std::cout << "read " << nLines << " lines" << std::endl;

  xSize = x_space.size();
  ySize = y_space.size();
  zSize = z_space.size();
}

void scalarField::set(int i, int j, int k, double value)
{
  values[ySize*zSize*i + zSize*j + k] = value;
}

double scalarField::get(int i, int j, int k)
{
  // std::cout << i << '\t' << j << '\t' << k << std::endl;
  return values[ySize*zSize*i + zSize*j + k];
}

void scalarField::print_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      for ( int k = 0; k < zSize; k++ ) {
	outFile << x_space[i] << ','
		<< y_space[j] << ','
		<< z_space[k] << ','
		<< get(i, j, k) << '\n';
      }
    }
  }
  outFile.close();
}

double scalarField::interpolate(std::vector <double> pos)
{
  std::vector <double> spacing = {0, 0, 0};
  
  // get the 8 nearest points on the grid
  int i = 0;
  while ( x_space[i] < pos[0] ) {
    i++;
  }
  int xLowInd = i-1;
  int xHighInd = i;

  spacing[0] = x_space[xHighInd] - x_space[xLowInd];
  
  i = 0;
  while ( y_space[i] < pos[1] ) {
    i++;
  }
  int yLowInd = i-1;
  int yHighInd = i;

  spacing[1] = x_space[yHighInd] - x_space[yLowInd];
  
  i = 0;
  while ( z_space[i] < pos[2] ) {
    i++;
  }
  int zLowInd = i-1;
  int zHighInd = i;

  spacing[2] = x_space[zHighInd] - x_space[zLowInd];
  
  // linear interpolation
  double sum = 0;
  double prod;
  for ( int i: {xLowInd, xHighInd} ) {
    for ( int j: {yLowInd, yHighInd} ) {
      for ( int k: {zLowInd, zHighInd} ) {
	std::vector <double> vertexPos = {x_space[i], y_space[j], z_space[k]};
	prod = get(i, j, k);
	for ( int n: {0, 1, 2} ) {
	  prod *= spacing[n] - abs( pos[n] - vertexPos[n] );
	}
	sum += prod;
      }
    }
  }

  for ( int n: {0, 1, 2} ) {
    sum /= spacing[n];
  }

  return sum;
  
  // inverse distance weighting
  // double normalization = 0;
  // double sum = 0;
  // for ( int i: {xLowInd, xHighInd} ) {
  //   for ( int j: {yLowInd, yHighInd} ) {
  //     for ( int k: {zLowInd, zHighInd} ) {
  // 	std::vector <double> edgePos = {x_space[i], y_space[j], z_space[k]};
  // 	double weight = pow(mag(pos - edgePos), -1);

  // 	sum += weight*get(i, j, k);
  // 	normalization += weight;
  // 	// std::cout << get(i, j, k) << std::endl;
  //     }
  //   }
  // }
  // sum /= normalization;
  
  // return sum;
}

std::vector <double> scalarField::interpolate_grad(std::vector <double> pos)
{
  std::vector <double> spacing = {0, 0, 0};

  std::vector <double> grad = {0, 0, 0};
  
  // get the 8 nearest points on the grid
  int i = 0;
  while ( x_space[i] < pos[0] ) {
    i++;
  }
  int xLowInd = i-1;
  int xHighInd = i;

  spacing[0] = x_space[xHighInd] - x_space[xLowInd];
  
  i = 0;
  while ( y_space[i] < pos[1] ) {
    i++;
  }
  int yLowInd = i-1;
  int yHighInd = i;

  spacing[1] = x_space[yHighInd] - x_space[yLowInd];
  
  i = 0;
  while ( z_space[i] < pos[2] ) {
    i++;
  }
  int zLowInd = i-1;
  int zHighInd = i;

  spacing[2] = x_space[zHighInd] - x_space[zLowInd];
  
  // x grad
  // interpolate in y and z
  double weight;
  double highProd;
  double lowProd;
  for ( int j: {yLowInd, yHighInd} ) {
    for ( int k: {zLowInd, zHighInd} ) {
      std::vector <double> vertexPos = {0, y_space[j], z_space[k]};
      highProd = get(xHighInd, j, k);
      lowProd = get(xLowInd, j, k);
      for ( int n: {1, 2} ) {
	weight = spacing[n] - abs(pos[n] - vertexPos[n]);
	highProd *= weight;
	lowProd *= weight;
      }
      grad[0] += highProd;
      grad[0] -= lowProd;
    }
  }

  // y grad
  // interpolate in x and z
  for ( int i: {xLowInd, xHighInd} ) {
    for ( int k: {zLowInd, zHighInd} ) {
      std::vector <double> vertexPos = {x_space[i], 0, z_space[k]};
      highProd = get(i, yHighInd, k);
      lowProd = get(i, yLowInd, k);
      for ( int n: {0, 2} ) {
	weight = spacing[n] - abs(pos[n] - vertexPos[n]);
	highProd *= weight;
	lowProd *= weight;
      }
      grad[1] += highProd;
      grad[1] -= lowProd;
    }
  }

  // z grad
  // interpolate in x and y
  for ( int i: {xLowInd, xHighInd} ) {
    for ( int j: {yLowInd, yHighInd} ) {
      std::vector <double> vertexPos = {x_space[i], y_space[j], 0};
      highProd = get(i, j, zHighInd);
      lowProd = get(i, j, zLowInd);
      for ( int n: {0, 1} ) {
	weight = spacing[n] - abs(pos[n] - vertexPos[n]);
	highProd *= weight;
	lowProd *= weight;
      }
      grad[2] += highProd;
      grad[2] -= lowProd;
    }
  }

  for ( int n: {0, 1, 2} ) {
    grad = grad*(1./spacing[n]);
  }

  return grad;
}

double squared_diff(scalarField a, scalarField b)
{
  // return the squared differences of each field

  // first, check that they have the same shape
  try {
    if ( ( not ( a.xSize == b.xSize ) )
	 or ( not ( a.ySize == b.ySize ) )
	 or ( not ( a.zSize == b.zSize ) ) ) {
      throw 20;
    }
    else {
      double sum = 0;
      int n_points = 0;
      for ( int i = 0; i < a.xSize; i++ ) {
	for ( int j = 0; j < a.ySize; j++ ) {
	  for ( int k = 0; k < a.zSize; k++ ) {
	    sum += pow(a.get(i, j, k) - b.get(i, j, k), 2);
	    n_points++;
	  }
	}
      }
      return sum/n_points;
    }
  }
  catch ( int e ) {
    std::cout << "ERROR: tried to compute squared difference of fields of different size!" << std::endl;
  }
  return 0;
}

double squared_diff(scalarField a, scalarField b, scalarField exclude)
{
  // return the squared differences of each field

  // first, check that they have the same shape
  try {
    if ( ( not ( a.xSize == b.xSize ) )
	 or ( not ( a.ySize == b.ySize ) )
	 or ( not ( a.zSize == b.zSize ) ) ) {
      throw 20;
    }
    else {
      double sum = 0;
      int n_points = 0;
      for ( int i = 0; i < a.xSize; i++ ) {
	for ( int j = 0; j < a.ySize; j++ ) {
	  for ( int k = 0; k < a.zSize; k++ ) {
	    if ( not exclude.get(i, j, k) ) {
	      sum += pow(a.get(i, j, k) - b.get(i, j, k), 2);
	      n_points++;
	    }
	  }
	}
      }
      return sum/n_points;
    }
  }
  catch ( int e ) {
    std::cout << "ERROR: tried to compute squared difference of fields of different size!" << std::endl;
  }
  return 0;
}
