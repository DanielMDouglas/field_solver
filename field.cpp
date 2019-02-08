#include <vector>
#include <string>
#include <fstream>

#include "field.h"

field::field(std::vector<double> x, std::vector<double> y, double init_val = 0) {
  xSize = x.size();
  ySize = y.size();

  x_space = x;
  y_space = y;

  values = std::vector< double > (xSize*ySize, init_val);
}

field::field(std::vector<double> x, std::vector<double> y, double (*f)(double, double)) {
  xSize = x.size();
  ySize = y.size();

  x_space = x;
  y_space = y;

  values = std::vector< double > (xSize*ySize);

  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      set(i, j, (*f)(x[i], y[j]));
    }
  }
}

void field::set(int i, int j, double value) {
  values[ySize*i + j] = value;
}

double field::get(int i, int j) {
  return values[ySize*i + j];
}

void field::print_to_file(std::string filename) {
  std::ofstream outFile (filename.c_str());
  for ( int i = 0; i < xSize; i++ ) {
    for ( int j = 0; j < ySize; j++ ) {
      outFile << x_space[i] << ','
	      << y_space[j] << ','
	      << get(i, j) << '\n';
    }
  }
  outFile.close();
}
