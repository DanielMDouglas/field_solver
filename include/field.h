#include "boundary.h"

class field
{
 public:
  int xSize, ySize, zSize;

  std::vector <double> values;
  std::vector <double> x_space;
  std::vector <double> y_space;
  std::vector <double> z_space;

  field(std::vector <double>,
	std::vector <double>,
	std::vector <double>,
	double);
  field(std::vector <double>,
	std::vector <double>,
	std::vector <double>,
	double (*)(double, double, double));
  field(boundary, int, int, int); 
  void set(int, int, int, double);
  double get(int, int, int);
  void print_to_file(std::string);
};

double squared_diff(field, field);
double squared_diff(field, field, bool (*)(int, int, int));
