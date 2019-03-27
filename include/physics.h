#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>

#include "scalar_field.h"
#include "path.h"

const double dt = 1.e-4; // us
const double Tb = 87.302; // boiling temperature of LAr
const double e = 1.603e-19; // C

std::vector <double> E(std::vector <double>, scalarField*);
std::vector <double> driftV(std::vector <double>, double);

/* path drift_path(std::vector <double>, scalarField*, boundary); */
void drift_path(std::vector <double>, scalarField*, boundary, path*&);
std::vector <double> ramo_induction(double, path *, scalarField *);
std::vector <double> sample(std::vector <double>);
void sample(std::vector <double>, double*, std::vector <double>*);

#endif
