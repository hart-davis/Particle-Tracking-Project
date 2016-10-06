import tracpy
import tracpy.plotting
#import mike_plotting
from tracpy.tracpy_class import Tracpy
import datetime
import matplotlib
from tracpy import op
#import Mike_Test_Case
#import tracpy_modified
#from tracpy import inout
import xray
from numpy import where, squeeze, array, sqrt, zeros
import os
import glob
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from numpy import empty, size, unique, diff, asfortranarray

matplotlib.rcParams.update({'font.size': 20})
# this is the location of a high resolution model TXLA (texas something) Change in all files
loc = '/home/hart-davis/Desktop/Mike_PTA/Globcurrent_data_final.nc'
# number of days that this model will be run. This example ranges from 1 to 24
ndays = 1
# this is the starting date in a datetime
date = datetime.datetime (2008, 3, 1, 0, 0, 0)
# this the time between outputs in the correct format (seconds)
tseas = 3*3600

time_units = 'seconds since 1970-01-01'
# sets a smaller limit than between model outputs for when to force interpolation if it hasn't already occured
nsteps = 400
# controls the sampling frequency of the drifter tracks
N = 400 #4
# controls the direction in time that the model is run 1 forward/ -1 backward
ff = 1

ah = 0
av = 0
# ranges from 0-3 and it changes the behaviour
doturb = 0
# this is just the name used to save the results into a netCDF file
name = 'temp'
# this is whether or not the tracker should be done in a 3D format (1-3D/0-2D)
do3d = 0

z0 = 's'

num_layers = 1

zpar = 0

# nwgom is a built i projection that is used. Also - galveston and nwgom-pyproj
proj = tracpy.tools.make_proj('madagascar')
# read in the grid
grid = tracpy.inout.readgrid(loc, proj,  usespherical = False)

grid.proj = proj
grid.sc_r = np.array([-1. ,0.])
grid.Cs_r = np.array([-1. ,0.])
grid.Vtransform = 1
grid.Vstretching = 1
grid.theta_s = 1
grid.theta_b = 0
grid.hc = 0
grid.km = 1

# initialization of the tracpy class
tp = Tracpy (loc, grid, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, time_units=time_units)
# input starting locations in the long and latitude
lon0, lat0 = np.meshgrid(np.linspace(43.125,47.5,30), np.linspace(-24.125,-25.5,30))
# removes points not found within the domain
lon0, lat0 = tracpy.tools.check_points(lon0, lat0, tp.grid)
# timing that the grid was already read in
lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)
# plot the tracks
tracpy.plotting.tracks(lonp, latp, tp.name, tp.grid)

#tracpy.plotting.hist(lonp, latp, tp.name, grid=tp.grid, which='pcolor', bins=(50,50))
