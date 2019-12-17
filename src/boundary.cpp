#include <iostream>
#include <fstream>
#include <vector>
#include <TRandom3.h>
#include <nlohmann/json.hpp>

#include "boundary.h"
#include "pad.h"

boundary::boundary(std::string inFileName)
{
  nlohmann::json j;
  std::ifstream i(inFileName);
  i >> j;

  std::cout << "# Building geometry: " << j["name"] << std::endl;

  periodicX = j["periodicity"]["x"];
  periodicY = j["periodicity"]["y"];
  periodicZ = j["periodicity"]["z"];
  
  Xmin = j["bounds"]["xmin"];
  Xmax = j["bounds"]["xmax"];
  Ymin = j["bounds"]["ymin"];
  Ymax = j["bounds"]["ymax"];
  Zmin = j["bounds"]["zmin"];
  Zmax = j["bounds"]["zmax"];

  for ( nlohmann::json volume_json : j["volumes"] ) {
    add_volume(new volume(volume_json));
  }
}

void boundary::add_volume(volume * newVol)
{
  volumes[nVolumes] = newVol;
  nVolumes++;
}

bool boundary::is_in_conductor(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> type == "conductor";
    }
  }
  return false;
}
  
bool boundary::is_in_VN(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> type == "von Neumann";
    }
  }
  return false;
}

bool boundary::is_in_volume(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return true;
    }
  }
  return false;
}

double boundary::boundary_value(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      if ( volumes[i] -> type == "conductor" ) {
	return volumes[i] -> V(x, y, z);
      }
      else if ( volumes[i] -> type == "von Neumann" ) {
	return 0xdeadbeef;
      }
    }
  }
  return 0;
}

double boundary::Efield(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( ( volumes[i] -> type == "von Neumann" )
	 and ( volumes[i] -> is_in_boundary(x, y, z) ) ) {
      return volumes[i] -> Efield;
    }
  }
  return 0;
}

double boundary::permittivity(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> er(x, y, z);
    }
  }
  return 1;
}

double boundary::conductivity(double x, double y, double z)
{
  for ( int i = 0; i < nVolumes; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      return volumes[i] -> sigma(x, y, z);
    }
  }
  return 0;
}
