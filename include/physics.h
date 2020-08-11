#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>
#include <math.h>

#include "field.h"
#include "path.h"

const double dt = 5.e-4; // us
const double Tb = 87.302; // boiling temperature of LAr
const double e = 1.603e-19; // C
const double ep0 = 8.85e-11; // C/(kV*cm)
const double pi = 4*atan(1);

std::vector <double> E(std::vector <double>, field <double>*);
std::vector <double> driftV(std::vector <double>, double);

/* path drift_path(std::vector <double>, field*, boundary); */
void drift_path(std::vector <double>, field <double>*, boundary, path*&);
std::vector <double> ramo_induction(double, path *, field <double>*);
std::vector <double> sample(std::vector <double>);
void sample(std::vector <double>, double*, std::vector <double>*);

#endif
