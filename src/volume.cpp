#include <iostream>

#include "volume.h"
#include "utils.h"

volume::volume(nlohmann::json volume_json)
{  
  std::cout << "# Constructing " << volume_json["type"]
	    << " volume, " << volume_json["description"] << std::endl;

  Xmin = volume_json["xmin"];
  Xmax = volume_json["xmax"];
  Ymin = volume_json["ymin"];
  Ymax = volume_json["ymax"];
  Zmin = volume_json["zmin"];
  Zmax = volume_json["zmax"];

  center = std::vector <double> {0.5*(Xmin + Xmax),
				 0.5*(Ymin + Ymax),
				 0.5*(Zmin + Zmax)};
  
  description = volume_json["description"];
  
  type = volume_json["type"];

  if ( type == "conductor" ) {
    er = constant(1);
    sigma = constant(0);
    if ( volume_json["voltage"]["function_name"] == "constant" ) {
      V = constant(volume_json["voltage"]["args"][0]);
    }
    else if ( volume_json["voltage"]["function_name"] == "linear" ) {
      V = linear(volume_json["voltage"]["args"][0],
		 volume_json["voltage"]["args"][1],
		 volume_json["voltage"]["args"][2],
		 volume_json["voltage"]["args"][3]);
    }
    else if ( volume_json["voltage"]["function_name"] == "gaussian" ) {
	V = gaussian(volume_json["voltage"]["args"][0],
		     volume_json["voltage"]["args"][1],
		     volume_json["voltage"]["args"][2],
		     volume_json["voltage"]["args"][3],
		     volume_json["voltage"]["args"][4],
		     volume_json["voltage"]["args"][5]);
    }
  }
  else if ( type == "von Neumann" ) {
    V = constant(0xdeadbeef);
    Efield = volume_json["Efield"];
    er = constant(1);
    sigma = constant(0);
  }
  else if ( type == "dielectric" ) {
    V = constant(0);
    if ( volume_json["permittivity"]["function_name"] == "constant" )
      {
	er = constant(volume_json["permittivity"]["args"][0]);
      }
    else if ( volume_json["permittivity"]["function_name"] == "linear" ) {
      er = linear(volume_json["permittivity"]["args"][0],
		  volume_json["permittivity"]["args"][1],
		  volume_json["permittivity"]["args"][2],
		  volume_json["permittivity"]["args"][3]);
    }
    else if ( volume_json["permittivity"]["function_name"] == "gaussian" ) {
      er = gaussian(volume_json["permittivity"]["args"][0],
		    volume_json["permittivity"]["args"][1],
		    volume_json["permittivity"]["args"][2],
		    volume_json["permittivity"]["args"][3],
		    volume_json["permittivity"]["args"][4],
		    volume_json["permittivity"]["args"][5]);
    }

    if ( volume_json["conductivity"]["function_name"] == "constant" ) {
      sigma = constant(volume_json["conductivity"]["args"][0]);
    }
    else if ( volume_json["conductivity"]["function_name"] == "linear" ) {
      sigma = linear(volume_json["conductivity"]["args"][0],
		     volume_json["conductivity"]["args"][1],
		     volume_json["conductivity"]["args"][2],
		     volume_json["conductivity"]["args"][3]);
    }
    else if ( volume_json["conductivity"]["function_name"] == "gaussian" ) {
      sigma = gaussian(volume_json["conductivity"]["args"][0],
		       volume_json["conductivity"]["args"][1],
		       volume_json["conductivity"]["args"][2],
		       volume_json["conductivity"]["args"][3],
		       volume_json["conductivity"]["args"][4],
		       volume_json["conductivity"]["args"][5]);
    }
  }
}

double volume::get_voltage(double x, double y, double z)
{
  return V(x, y, z);
}

bool volume::is_in_boundary(double x, double y, double z)
{
  if ( ( Xmin <= x ) and
       ( x <= Xmax ) and
       ( Ymin <= y ) and
       ( y <= Ymax ) and
       ( Zmin <= z ) and
       ( z <= Zmax ) ) {
    return true;
  }
  else return false;
}
