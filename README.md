# Electrostatic Field Solver for LArTPC's

## Overview

This package contains software for simulating static electric fields and the properties of electrons drifting within them

## Dependencies

This package relies on nlohmann's JSON package (https://github.com/nlohmann/json)

## Installation

Installation is handled through `cmake`:

```
mkdir build
cd build
cmake ..
make
make install
```

## Executables

### solver

This program uses the finite difference method (FDM) to solve the electrostatic potential inside of a volume.  For more details about the mathematical methods, please see James R. Nagel's "Solving the Generalized Poisson Equation Using the Finite-Difference Method (FDM)".

```
Usage: ./solver [OPTIONS]
-i      starting solution, default: none
-q      starting charge distribution, default: none
-o      output, default: final.dat
-g      geometry JSON, default: geometries/bulkPix.json
-n      maximum number of iterations, default: 100000
-f      frequency of reporting, default: 100
-t      threshold for stopping relaxation, default: 1e-15
-w      over-relaxation factor, default: 1
-N      number of vertices, default: number specified by spacing
-u      number of upscaling operations, default: 0
-s      grid spacing, default: 0.01 (cm)
-v      verbosity (each reporting on a new line), default: false
-h      display this help and exit
```

## Getting Started

As a very basic example, let's run through a drift inside of a very simple geometry.  First, we need to solve the potential field:

```
./solver -g geometries/linear.json -s 0.2
```

This will use the geometry defined in `linear.json` and solve it on a grid with spacing of 0.2 units

## Geometries

Let's look at the makeup of a geometry description.  As above, we'll use `linear.json`, as it's pretty simple:

```
{
    "name": "linear",
    "periodicity": {"x": true,
	            "y": true,
		    "z": false
		    },
    "bounds": {"xmin": 0,
    	       "xmax": 1,
	       "ymin": 0,
	       "ymax": 1,
	       "zmin": 0,
	       "zmax": 1
	       },
    "volumes": [{"description": "bottom plate",
		 "type": "conductor",
		 "xmin": 0,
		 "xmax": 1,
		 "ymin": 0,
		 "ymax": 1,
		 "zmin": 0,
		 "zmax": 0.01,
		 "voltage": {"function_name": "constant",
		 	     "args": [0]
		            }
		},
		{"description": "top plate",
		 "type": "conductor",
		 "xmin": 0,
		 "xmax": 1,
		 "ymin": 0,
		 "ymax": 1,
		 "zmin": 0.99,
		 "zmax": 1,
		 "voltage": {"function_name": "constant",
		             "args": [1]
			    }
		},
		{"description": "bulk volume",
		 "type": "dielectric",
		 "xmin": 0,
		 "xmax": 1,
		 "ymin": 0,
		 "ymax": 1,
		 "zmin": 0,
		 "zmax": 1,
		 "permittivity": {"function_name": "constant",
					           "args": [1]
			         },
		 "conductivity": {"function_name": "constant",
			 	  "args": [1]
				  }
                 }
		]
}
```