#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

def read_csv(infileName):
    data = np.loadtxt(infileName, delimiter=',')
    x = np.unique(data[:,0])
    y = np.unique(data[:,1])
    z = np.unique(data[:,2])
    values = np.reshape(data[:,3], (x.size, y.size, z.size), order='A').T
    
    return x, y, z, values

def choose_slice(x, y, z, values, choice):
    """
    pass the data as well as a 'choice' string
    the string should be similar to 
    z = 5
    or
    x = 2
    which defines a plane in units of pixel number
    """

    choice_var = choice.split('=')[0].strip()
    choice_val = int(choice.split('=')[1].strip())

    if choice_var == 'x':
        u = y
        v = z
        selec = values[choice_val,:,:]
    elif choice_var == 'y':
        u = z
        v = x
        selec = values[:,choice_val,:]
    elif choice_var == 'z':
        u = x
        v = y
        selec = values[:,:,choice_val]

    return u, v, selec

if __name__ == '__main__':
    infileNames = sys.argv[1:-1]
    choice_string = sys.argv[-1]

    if not infileNames:
        print "Please specify the file to plot from!"
        print "Usage: plot [fieldData.dat]"

        exit

    for infileName in infileNames:
        x, y, z, values = read_csv(infileName)
        u, v, selec = choose_slice(x, y, z, values, choice_string)

        Vmax = np.max(selec)
    
        plt.figure()
        cmesh = plt.pcolormesh(u, v, selec, cmap = 'Blues')
        cb = plt.colorbar(ticks = [0, Vmax])
        cmesh.set_clim(vmin=0, vmax=Vmax)
        cb.ax.set_yticklabels([0, r'$V_{max}$'])

        plt.contour(u, v, selec, 20)
        
        plt.xlabel(r'$x$ [m]')
        plt.ylabel(r'$y$ [m]')

        plt.show()
