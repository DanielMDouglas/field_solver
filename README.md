# Electrostatic Field Solver for LArTPC's

## Overview

This package contains software for simulating static electric fields and the properties of electrons drifting within them

## Dependencies

This package relies on nlohmann's JSON package (https://github.com/nlohmann/json).

To use many of the visualization tools in the `utils` directory, you will also need Matplotlib and Numpy.

## Installation

Installation is handled through `cmake`:

```
mkdir build
cd build
cmake ..
make
make install
```

## Getting Started

As a very basic example, let's run through a drift inside of a very simple geometry.  First, we need to solve the potential field.  From the `build` directory:

```
./solver -g geometries/linear.json -s 0.2
```

This will use the geometry defined in `linear.json` and solve it on a grid with spacing of 0.2 units.  The result is saved to `final.dat` by default.

If we want to take a look at the potential we just solved, we can do

```
./utils/plot final.dat 'y = 3' 'Potential'
```

This should show us our coarse grid over the plane y = 3 and give it the title 'Potential' (TODO: add the image to this README).

Next, we need to calculate the weighting potential.  This is the potential as it would be if the electrod we are interested in (here it's the top plate) was at a unit potential and everything else is set to ground.  Fortunately, `linear.json` is already configured this way, so in this case the potential and the weighting potential are identical!  In other cases, it is usually necessary to define a separate geometry JSON with the potentials appropriately tweaked.

The next step is to calculate the path that a drifting charge will take within this field:

```
./drift_single -x 0.5 -y 0.5 -z 0.5 -o drift.dat -f final.dat -g geometries/linear.json
```

This will place an electron at (x, y, z) = (0.5, 0.5, 0.5) and let it drift until its position is within a conductor or otherwise outside of the defined volume.  The output containing time, position, and velocity information is saved to `drift.dat`.

The last step is using this drift path and the weighting potential to calculate the current induced in the top plate as the charge drifts towards it:

```
./induction_single -d drift.dat -w final.dat -o current.dat
```

This step uses the Ramo Theorem (https://en.wikipedia.org/wiki/Shockley%E2%80%93Ramo_theorem).  The resulting output, saved to `current.dat`, is a time-series of current.

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

### drift_single

This program allows one to drift single charges at a time from a given position.

```
Usage: ./drift_single [OPTIONS]
-x      Initial x position, default: 0
-y      Initial y position, default: 0
-z      Initial z position, default: 0.45
-o      output, default: drift.dat
-f      electric potential field, default: none
-g      detector geometry, default: geometries/bulkPix.json
-h      display this help and exit
```

### induction_single

This program uses a pre-calculated path and a weighting field to produce a current series using the Ramo Theorem.

```
Usage: ./induction_single [OPTIONS]
-d      potential field, default: none
-w      weighting field, default: none
-o      output, default: none
-h      display this help and exit
```

## Geometries

Let's look at the makeup of a geometry description.  As above, we'll use `linear.json`, as it's pretty simple:

```
{
    "name": "linear",                                     # a simple descriptor, can be anything
    "periodicity": {"x": true,			          # which bounds are periodic, if any
	            "y": true,
		    "z": false
		    },
    "bounds": {"xmin": 0,			          # define the extent of the overall geometry
    	       "xmax": 1,
	       "ymin": 0,
	       "ymax": 1,
	       "zmin": 0,
	       "zmax": 1
	       },
    "volumes": [{"description": "bottom plate",           # each volume has a name
		 "type": "conductor",		          # volume type
		 "xmin": 0,
		 "xmax": 1,
		 "ymin": 0,
		 "ymax": 1,
		 "zmin": 0,
		 "zmax": 0.01,
		 "voltage": {"function_name": "constant", # for conductors, you must define the voltage
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

### Volume types

volumes can be `"conductor"`, `"dielectric"` or `"neumann"`.

`"conductor"` volumes are conductors, essentially dielectrics with infinite conductivity.  These volumes must also contain a "voltage" key with a functional description of the voltage.

`"dielectric"` volumes are normal matter, typically liquid Argon, plastic, or something similar.  These volumes must contain a "permittivity" and a "conductivity" key, each with functional descriptions.

`"neumann"` volumes are not quite real, physical volumes, but they are useful for fixing the electric field normal to a boundary. These volumes must contain an "Efield" key and a simple float number to define the normal component of the electric field.

### Functional descriptors

functions (as of the writing of this README) can be `"constant"`, `"linear"`, or `"gaussian"`.  Arguments are passed as a list with the key `"args"`.

#### `"constant"`

| Argument Position | Meaning |
| ----------------- | ------- |
| 0                 | value   |

#### `"linear"`

| Argument Position | Meaning   |
| ----------------- | --------- |
| 0                 | intercept |
| 1                 | x slope   |
| 2                 | y slope   |
| 3                 | z slope   |

#### `"gaussian"`

| Argument Position | Meaning    |
| ----------------- | ---------- |
| 0                 | x center   |
| 1                 | y center   |
| 2                 | z center   |
| 3                 | width      |
| 4                 | baseline   |
| 5                 | magnitude  |
