#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <math.h>
#include <boost/algorithm/string.hpp>

#include <TFile.h>
#include <TH3.h>

// #include "field.h"

template <typename T>
field<T>::field(std::vector <double> x,
		std::vector <double> y,
		std::vector <double> z,
		T init_val)
{
  xSize = x.size();
  ySize = y.size();
  zSize = z.size();

  x_space = x;
  y_space = y;
  z_space = z;

  // values = std::vector< double > (xSize*ySize*zSize, init_val);
  // values = std::array <double, size>;

  values = new T [xSize*ySize*zSize];
  
  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      for ( int k = 0; k < zSize; k++ ) {
	set(i, j, k, init_val);
      }
    }
  }
}

template <typename T>
field<T>::field(std::vector <double> x,
		std::vector <double> y,
		std::vector <double> z,
		std::function <T (double, double, double)> f)
{
  xSize = x.size();
  ySize = y.size();
  zSize = z.size();

  x_space = x;
  y_space = y;
  z_space = z;

  // values = new std::array <T, xSize*ySize*zSize>;
  // allocate(xSize*ySize*zSize);
  values = new T [xSize*ySize*zSize];
  
  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      for ( int k = 0; k < zSize; k++ ) {
	set(i, j, k, f(x[i], y[j], z[k]));
      }
    }
  }
}

template <typename T>
field<T>::field(std::string filename)
{
  std::vector<std::string> split;
 
  boost::split(split, filename, [](char c){return c == '.';});

  std::cout << "# Reading from " << filename << std::endl;
  // std::cout << "which has extension " << split[1] << std::endl;
  if ( split[1] == "root" ) {
    // load root file
    TFile * inFile = new TFile(filename.c_str());

    // TH3D * hist = dynamic_cast <TH3D*> (inFile -> Get("field"));
    TH3D * hist = static_cast <TH3D*> (inFile -> Get("field"));
    
    xSize = hist -> GetNbinsX();
    x_space = std::vector <double> ();
    for ( int i = 0; i < xSize; i++ ) {
      x_space.push_back(hist -> GetXaxis() -> GetBinCenter(i));
    }
    
    ySize = hist -> GetNbinsY();
    y_space = std::vector <double> ();
    for ( int j = 0; j < ySize; j++ ) {
      y_space.push_back(hist -> GetYaxis() -> GetBinCenter(j));
    }
    
    zSize = hist -> GetNbinsZ();
    z_space = std::vector <double> ();
    for ( int k = 0; k < ySize; k++ ) {
      z_space.push_back(hist -> GetZaxis() -> GetBinCenter(k));
    }

    values = new T [xSize*ySize*zSize];
    
    for ( int i = 0; i < xSize; i++ ) {
      for ( int j = 0; j < ySize; j++ ) {
	for ( int k = 0; k < zSize; k++ ) {
	  set(i, j, k, hist -> GetBinContent(i, j, k));
	}
      }
    }

    inFile -> Close();
  }
  else {
    // treat it as a plaintext file
    std::ifstream inFile (filename.c_str());

    std::string entry;
    double value;
    int nLines = 0;

    std::vector <T> vals;

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
      vals.push_back(std::stod(entry));
    
      nLines++;
    }
  
    inFile.close();
    std::cout << "# read " << nLines << " lines" << std::endl;
  
    xSize = x_space.size();
    ySize = y_space.size();
    zSize = z_space.size();

    // values = new std::array <T, nLines>;
    values = new T [nLines];
  
    for ( int line_iter = 0; line_iter < nLines; line_iter++ ) {
      values[line_iter] = vals[line_iter];
    }
  }
}

template <typename T>
void field<T>::set(int i, int j, int k, T value)
{
  values[ySize*zSize*i + zSize*j + k] = value;
}

template <typename T>
T field<T>::get(int i, int j, int k)
{
  if ( ( i == 0xdeadbeef ) or ( j == 0xdeadbeef ) or ( k == 0xdeadbeef ) ) {
    // if ( T == std::vector ) {
    //   return {0, 0, 0};
    // }
    // else {
    //   return 0;
    // }
    return T (0);
  }
  else {
    return values[ySize*zSize*i + zSize*j + k];
  }
}

template <typename T>
void field<T>::print_to_file(std::string filename)
{
  std::vector <std::string> split;
  boost::split(split, filename, [](char c){return c == '.';});
  std::cout << "# Writing values to " << filename << std::endl;

  if ( split[1] == "root" ) {
    // write to root file
    TFile * outFile = new TFile(filename.c_str(), "NEW");
    
    TH3D * hist = new TH3D ("field", "field",
			    xSize, x_space[0], x_space[xSize-1],
			    ySize, y_space[0], y_space[ySize-1],
			    zSize, z_space[0], z_space[zSize-1]);
    for ( int i = 0; i < xSize; i++ ) {
      for ( int j = 0; j < ySize; j++ ) {
	for ( int k = 0; k < zSize; k++ ) {
	  hist -> SetBinContent(i, j, k, get(i, j, k));
	}
      }
    }

    hist -> Write();
    outFile -> Close();
  }
  else {
    // treat it as a plaintext file
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
}

template <typename T>
field<T> * field<T>::upscale(int scaleFactor)
{
  // make the grid finer by scaleFactor
  // keep the same endpoints
  // so the number of points isn't just e.g. xSize*scaleFactor
  // really, it's dividing the spacing
  float xMin = x_space[0];
  float yMin = y_space[0];
  float zMin = z_space[0];

  float xMax = x_space[xSize - 1];
  float yMax = y_space[ySize - 1];
  float zMax = z_space[zSize - 1];
    
  int newXSize = scaleFactor*(xSize - 1) + 1;
  int newYSize = scaleFactor*(ySize - 1) + 1;
  int newZSize = scaleFactor*(zSize - 1) + 1;

  field<T> * newField = new field<T> (linspace(xMin, xMax, newXSize),
				      linspace(yMin, yMax, newYSize),
				      linspace(zMin, zMax, newZSize),
				      0);

  int sub_i_endpoint;
  int sub_j_endpoint;
  int sub_k_endpoint;
  
  for ( int i = 0; i < xSize; i++ ) {
    if ( i == xSize - 1 ) {
      sub_i_endpoint = 1;
    }
    else {
      sub_i_endpoint = scaleFactor;
    }
    for ( int sub_i = 0; sub_i < sub_i_endpoint; sub_i++ ) {
      for ( int j = 0; j < ySize; j++ ) {
	if ( j == ySize - 1 ) {
	  sub_j_endpoint = 1;
	}
	else {
	  sub_j_endpoint = scaleFactor;
	}
	for ( int sub_j = 0; sub_i < sub_j_endpoint; sub_j++ ) {
	  for ( int k = 0; k < zSize; k++ ) {
	    if ( k == zSize - 1 ) {
	      sub_k_endpoint = 1;
	    }
	    else {
	      sub_k_endpoint = scaleFactor;
	    }
	    for ( int sub_k = 0; sub_k < sub_k_endpoint; sub_k++ ) {
	      newField -> set(scaleFactor*i + sub_i,
			      scaleFactor*j + sub_j,
			      scaleFactor*k + sub_k,
			      get(i, j, k));
	    }
	  }
	}
      }
    }
  }

  return newField;
}
	
template <typename T>
T field<T>::interpolate(std::vector <double> pos)
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
  T sum = 0;
  T prod;
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
}

template <typename T>
std::vector <T> field<T>::interpolate_grad(std::vector <double> pos)
{
  std::vector <double> spacing = {0, 0, 0};

  std::vector <T> grad = {0, 0, 0};
  
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

  spacing[1] = y_space[yHighInd] - y_space[yLowInd];
  
  i = 0;
  while ( z_space[i] < pos[2] ) {
    i++;
  }
  int zLowInd = i-1;
  int zHighInd = i;

  spacing[2] = z_space[zHighInd] - z_space[zLowInd];
  
  // x grad
  // interpolate in y and z
  T weight;
  T highProd;
  T lowProd;
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
