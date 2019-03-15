#ifndef DRIFT_H
#define DRIFT_H

#include <vector>

#include "scalar_field.h"
#include "path.h"
#include "boundary.h"
#include "utils.h"
#include "physics.h"

void drift_and_save(double, scalarField*, boundary);

#endif
