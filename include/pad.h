#ifndef PAD_H
#define PAD_H

#include <vector>
#include <string>

#include "field.h"
#include "path.h"

class pad
{
 public:
  double total_time;
  std::vector <double> center;
  std::vector <double> current;
  std::vector <std::vector <double>> output;

  pad(double, std::vector <double>);
  void add_response(double,
		    double,
		    std::vector <double>,
		    path *,
		    field <double> *,
		    field <double> *,
		    boundary);
  void calculate_output();
  void print_current_to_file(std::string);
  void print_output_to_file(std::string);
};

#endif
