#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <vector>
#include <sstream>

using arg_t = std::stringstream *;
using fn_t = std::function <void (arg_t)>;

class argParser
{
 public:
  std::vector <std::pair <std::string, fn_t>> options;
  std::vector <std::pair <std::string, fn_t>> flags;

  void add_option(std::string, fn_t);
  void add_flag(std::string, fn_t);
  void parse(int, const char **);
};

#endif
