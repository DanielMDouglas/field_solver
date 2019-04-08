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
  std::vector <double> response;

  pad(double, std::vector <double>);
  void add_response(double,
		    double,
		    std::vector <double>,
		    field <double> *,
		    field <double> *,
		    boundary);
  void add_response(double,
		    double,
		    std::vector <double>,
		    path *,
		    field <double> *,
		    field <double> *,
		    boundary);
  void print_to_file(std::string);
};

#endif
