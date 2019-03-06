#ifndef DRIFT_H
#define DRIFT_H

#include "scalar_field.h"
#include "path.h"
#include "boundary.h"
#include "utils.h"

path drift_path(std::vector <double>, scalarField, boundary);

#endif
