#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include <string>

#include "volume.h"
#include "utils.h"

class boundary
{
 public:
  const static int maxNVolumes = 200;
  int nVolumes = 0;
  volume * volumes[maxNVolumes];
  double Xmin, Xmax = 0;
  double Ymin, Ymax = 0;
  double Zmin, Zmax = 0;
  bool periodicX = false;
  bool periodicY = false;
  bool periodicZ = false;
  
  boundary(std::string);
  void make_linear();
  void make_linear_cond_defect();
  void make_linear_diel_defect();
  void make_sheet();
  void make_sheet_cond_defect();
  void make_sheet_random_cond_defect();
  void make_sheet_diel_defect();
  void make_capacitor();
  void make_capacitor_with_dielectric();
  void make_bulkPix();
  void make_bulkPix_single();
  void make_bulkPixWeighting();
  void make_bulkWires();
  void add_volume(volume*);
  bool is_in_boundary(double, double, double);
  bool is_in_volume(double, double, double);
  double boundary_value(double, double, double);
  double permittivity(double, double, double);
  double conductivity(double, double, double);
};

#endif
