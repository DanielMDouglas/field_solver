#ifndef PAD_H
#define PAD_H

#include <vector>
#include <string>

#include "scalar_field.h"
#include "path.h"

class pad
{
 public:
  double total_time;
  std::vector <double> center;
  std::vector <double> response;

  pad(double, std::vector <double>);
  void add_response(double,
		    double,
		    std::vector <double>,
		    scalarField <double> *,
		    scalarField <double> *,
		    boundary);
  void add_response(double,
		    double,
		    std::vector <double>,
		    path *,
		    scalarField <double> *,
		    scalarField <double> *,
		    boundary);
  void print_to_file(std::string);
};

#endif
