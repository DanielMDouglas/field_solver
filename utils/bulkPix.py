import numpy as np
from json import encoder

wall_thickness = 0.02
pad_pitch = 0.4
pad_width = 0.3

xmin = 0 - wall_thickness
ymin = 0 - wall_thickness

zmin = 0
zmax = 0.5

def generate_central_pad(voltage):
    return {"type": "conductor",
            "xmin": xmin,
            "xmax": pad_width/2,
            "ymin": ymin,
            "ymax": pad_width/2,
            "zmin": zmin,
            "zmax": zmin + wall_thickness,
            "voltage": {"function_name": "constant",
                        "args": [voltage]},
            "description": "central pad"}

def generate_boundaries(xmax, ymax, weighting = False):
    if weighting:
        voltage = 0.
        Efield = 0.
    else:
        voltage = -20.
        Efield = 500.
    boundaries = [{"type": "conductor",
                   "xmin": xmin,
                   "xmax": xmax,
                   "ymin": ymin,
                   "ymax": ymax,
                   "zmin": zmin,
                   "zmax": zmin + wall_thickness,
                   "voltage": {"function_name": "constant",
                               "args": [voltage]},
                   "description": "back plane"},
                  {"type": "von Neumann",
                   "xmin": xmin,
                   "xmax": xmax,
                   "ymin": ymin,
                   "ymax": ymax,
                   "zmin": zmax - wall_thickness,
                   "zmax": zmax,
                   "Efield": Efield,
                   "description": "to main volume"},
                  {"type": "von Neumann",
                   "xmin": xmin,
                   "xmax": xmin + wall_thickness,
                   "ymin": ymin,
                   "ymax": ymax,
                   "zmin": zmin,
                   "zmax": zmax,
                   "Efield": 0,
                   "description": "x symmetry"},
                  {"type": "von Neumann",
                   "xmin": xmin,
                   "xmax": xmax,
                   "ymin": ymin,
                   "ymax": ymin + wall_thickness,
                   "zmin": zmin,
                   "zmax": zmax,
                   "Efield": 0,
                   "description": "y symmetry"},
                  {"type": "von Neumann",
                   "xmin": xmax - wall_thickness,
                   "xmax": xmax,
                   "ymin": ymin,
                   "ymax": ymax,
                   "zmin": zmin,
                   "zmax": zmax,
                   "Efield": 0,
                   "description": "x symmetry"},
                  {"type": "von Neumann",
                   "xmin": xmin,
                   "xmax": xmax,
                   "ymin": ymax - wall_thickness,
                   "ymax": ymax,
                   "zmin": zmin,
                   "zmax": zmax,
                   "Efield": 0,
                   "description": "y symmetry"}]
    return boundaries

def generate_other_pads(nPads):
    pads = []
    for xi in range(nPads):
        for yi in range(nPads):
            if not xi == yi == 0:
                x_center = pad_pitch*xi
                y_center = pad_pitch*yi
                
                pads.append({"type": "conductor",
                             "xmin": x_center - pad_width/2,
                             "xmax": x_center + pad_width/2,
                             "ymin": y_center - pad_width/2,
                             "ymax": y_center + pad_width/2,
                             "zmin": zmin,
                             "zmax": zmin + wall_thickness,
                             "voltage": {"function_name": "constant",
                                         "args": [0.]},
                             "description": "other pad"})

    return pads

def generate_argon(xmax, ymax):
    return {"type": "dielectric",
            "xmin": xmin,
            "xmax": xmax,
            "ymin": ymin,
            "ymax": ymax,
            "zmin": zmin,
            "zmax": zmax,
            "permittivity": {"function_name": "constant",
                             "args": [4]},
            "conductivity": {"function_name": "constant",
                             "args": [0]},
            "description": "argon"}
