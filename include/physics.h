#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>

#include "scalar_field.h"

std::vector <double> E(std::vector <double>, scalarField*);
std::vector <double> driftV(std::vector <double>, double);

#endif
