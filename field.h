/* #include <vector> */
/* #include <string> */

class field {
 public:
  int xSize, ySize;

  std::vector <double> values;
  std::vector <double> x_space;
  std::vector <double> y_space;

  field(std::vector <double>, std::vector <double>, double);
  field(std::vector <double>, std::vector <double>, double (*)(double, double));
  void set(int, int, double);
  double get(int, int);
  void print_to_file(std::string);
};
