#include <iostream>
#include <sstream>
#include <vector>
#include <functional>

#include "argParser.h"

void argParser::add_option(std::string opt, fn_t func)
{
  options.push_back(std::pair<std::string, fn_t> (opt, func));
}

void argParser::add_flag(std::string flag, fn_t func)
{
  flags.push_back(std::pair<std::string, fn_t> (flag, func));
}

void argParser::parse(int argc, const char ** argv)
{
  int opt = 0;
  while ( opt < argc ) {
    std::stringstream optValue;
    std::stringstream argValue;
    optValue << argv[++opt];
    argValue << argv[++opt];

    for ( std::pair <std::string, fn_t> thisOption : options ) {
      if ( optValue.str() == thisOption.first ) {
	thisOption.second(& argValue);
      }
    }
    for ( std::pair <std::string, fn_t> thisOption : options ) {
      if ( optValue.str() == thisOption.first ) {
	thisOption.second(& argValue);
	opt--;
      }
    }
  }
}
