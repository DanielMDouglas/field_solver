#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "field.h"
#include "solver.h"

field::field(std::vector <double> x,
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

field::field(std::vector <double> x,
	     std::vector <double> y,
	     std::vector <double> z,
	     double (*f)(double, double, double))
{
  xSize = x.size();
  ySize = y.size();
  zSize = z.size();

  x_space = x;
  y_space = y;
  z_space = z;

  values = std::vector< double > (xSize*ySize*zSize);

  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      for ( int k = 0; k < zSize; k++ ) {
	set(i, j, k, (*f)(x[i], y[j], z[k]));
      }
    }
  }
}

field::field(boundary bound, int Nx, int Ny, int Nz)
{
  xSize = Nx;
  ySize = Ny;
  zSize = Nz;

  x_space = linspace(bound.Xmin, bound.Xmax, Nx);
  y_space = linspace(bound.Ymin, bound.Ymax, Ny);
  z_space = linspace(bound.Zmin, bound.Zmax, Nz);

  values = std::vector< double > (xSize*ySize*zSize);

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

void field::set(int i, int j, int k, double value)
{
  values[ySize*zSize*i + zSize*j + k] = value;
}

double field::get(int i, int j, int k)
{
  return values[ySize*zSize*i + zSize*j + k];
}

void field::print_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      for ( int k = 0; k < zSize; k++ ) {
	outFile << x_space[i] << ','
		<< y_space[j] << ','
		<< z_space[k] << ','
		<< get(i, j, k) << '\n';
	// std::cout << x_space[i] << ','
	// 	  << y_space[j] << ','
	// 	  << z_space[k] << ','
	// 	  << get(i, j, k) << std::endl;;
      }
    }
  }
  outFile.close();
}

double squared_diff(field a, field b)
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

double squared_diff(field a, field b, bool (*exclude)(int, int, int))
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
	    if ( not exclude(i, j, k) ) {
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
