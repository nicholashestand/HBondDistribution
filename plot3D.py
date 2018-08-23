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
    parser.add_argument('files', metavar='File', nargs='+', 
                    help='the name of the file(s) containing the data')
    return parser
##################################################
       
# get the user input
args = get_user_args()

# store the user arguments in variables
files     = args.parse_args().files[0]
data      = np.loadtxt(files)
dimx      = len(np.unique(data[:,0]))
dimy      = len(np.unique(data[:,1]))
x = np.array(data[:,0]).reshape(dimx,dimy)
y = np.array(data[:,1]).reshape(dimx,dimy)
z = np.array(data[:,2]).reshape(dimx,dimy)


fig, ax = plt.subplots(1, figsize=(4,3), dpi=300 )
CS = ax.contourf( x, y, z, 200, cmap=plt.cm.jet )

ax.set_xlabel( r"R (nm)", fontsize=8)
ax.set_ylabel( r"$\beta$ (deg)", fontsize=8 )
ax.tick_params(axis='both',which='major',labelsize=8, direction='in',
                top='on',bottom='on', left='on', right='on')
ax.minorticks_on()
ax.tick_params(axis='both',which='minor',labelsize=6, direction='in',
                top='on',bottom='on', left='on', right='on')
ax.xaxis.set_major_locator(MultipleLocator(.1))
ax.xaxis.set_minor_locator(MultipleLocator(.02))
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(MultipleLocator(4))
ax.set_xlim((0.2,0.6))
ax.set_ylim((0,100))
fig.colorbar(CS,ticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])


plt.tight_layout()
plt.savefig("plt.eps")
