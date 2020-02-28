#ifndef DRIFT_H
#define DRIFT_H

#include <vector>

#include "field.h"
#include "path.h"
#include "boundary.h"
#include "utils.h"
#include "physics.h"

/* void drift_and_save(double, field*, boundary); */
void handleOpts(int, const char **);
void saySettings();
void drift_and_save(double, field <double>*, boundary, path*);

#endif
