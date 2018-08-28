#!/usr/bin/env python

# Read in a matrix of data a plot a specified line
import matplotlib.pyplot as plt
import numpy as np
import argparse 
from matplotlib.ticker import MultipleLocator

def get_user_args():
    """Read command line arguments and return as an argparse type""" 
    parser = argparse.ArgumentParser(description='Plot file(s) specified by the user. '+\
             'The format of the files should be space delimated columns containing the x, '+\
             'y and yerror data. Note that all of these options may not be available for zplot yet!')
    parser.add_argument('files', metavar='File', nargs=1, 
                    help='the name of the file(s) containing the data')
    parser.add_argument('--output', '-o', metavar='File', nargs=1, 
                    help='the name of the output file')
    
    return parser
##################################################
       
# get the user input
args = get_user_args()

# store the user arguments in variables
files     = args.parse_args().files[0]
output    = args.parse_args().output[0]
data      = np.loadtxt(files)
dimx      = len(np.unique(data[:,0]))
dimy      = len(np.unique(data[:,1]))
x = np.array(data[:,0]).reshape(dimx,dimy)
y = np.array(data[:,1]).reshape(dimx,dimy)
z = np.array(data[:,2]).reshape(dimx,dimy)

# saddle point at 300 K appears to be at 37 deg, 0.32
for i in np.argwhere( x==0.321 ):
    for j in np.argwhere( y==37.25 ):
        if np.all(i == j):
            saddlePoint = z[i[0],i[1]]

# normalize -- set saddle point to zero
z -= saddlePoint

# cut off the high values -- not sure how else to fix the colorbar here
for i in range(dimx):
    for j in range(dimy):
        if ( z[i,j] >= 3 and z[i,j] !=np.inf ): z[i,j] = 3

# make the plot
fig, ax = plt.subplots(1, figsize=(4,3), dpi=300 )
CS = ax.contourf( x, y, z, 24, cmap=plt.cm.jet, vmin=-3, vmax=3, extend='both' )
# get rid of the contour lines
#for c in CS.collections:
    #c.set_edgecolor("face")

# set the axis labels
ax.set_xlabel( r"R (nanometers)", fontsize=10)
ax.set_ylabel( r"$\beta$ (degrees)", fontsize=10 )
ax.tick_params(axis='both',which='major',labelsize=10, direction='in',
                top='on',bottom='on', left='on', right='on')
ax.minorticks_on()
ax.tick_params(axis='both',which='minor',labelsize=8, direction='in',
                top='on',bottom='on', left='on', right='on')
ax.xaxis.set_major_locator(MultipleLocator(.1))
ax.xaxis.set_minor_locator(MultipleLocator(.02))
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(MultipleLocator(4))
ax.set_xlim((0.2,0.6))
ax.set_ylim((0,100))
fig.colorbar(CS,ticks=MultipleLocator(1), extend='both')


plt.tight_layout()
plt.savefig(output)
plt.show()
