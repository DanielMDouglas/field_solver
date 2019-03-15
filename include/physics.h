#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>

#include "scalar_field.h"
#include "path.h"

std::vector <double> E(std::vector <double>, scalarField*);
std::vector <double> driftV(std::vector <double>, double);

path drift_path(std::vector <double>, scalarField*, boundary);
std::vector <double> ramo_induction(path, scalarField *);

#endif
