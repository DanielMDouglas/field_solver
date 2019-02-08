#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

def read_csv(infileName):
    data = np.loadtxt(infileName, delimiter=',')
    x = np.unique(data[:,0])
    y = np.unique(data[:,1])
    values = np.reshape(data[:,2], (x.size, y.size), order='A').T
    # x = data[:,0]
    # y = data[:,1]
    # values = data[:,2]
    
    return x, y, values
    
if __name__ == '__main__':
    infileNames = sys.argv[1:]

    if not infileNames:
        print "Please specify the file to plot from!"
        print "Usage: plot [fieldData.dat]"

        exit

    for infileName in infileNames:
        x, y, values = read_csv(infileName)

        Vmax = np.max(values)
    
        plt.figure()
        cmesh = plt.pcolormesh(x, y, values, cmap = 'Blues')
        cb = plt.colorbar(ticks = [0, Vmax])
        cmesh.set_clim(vmin=0, vmax=Vmax)
        cb.ax.set_yticklabels([0, r'$V_{max}$'])

        plt.contour(x, y, values, 20)
        
        plt.xlabel(r'$x$ [m]')
        plt.ylabel(r'$y$ [m]')

        plt.show()
