# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:54:23 2016

@author: hart-davis
"""

#I used iPython
#modules that are required for the code to run
import numpy as np
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
from tracpy.tracpy_class import Tracpy
import datetime
import matplotlib

matplotlib.rcParams.update({'font.size': 10})
# this is the location of a high resolution model TXLA (texas something) Change in all files
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# number of days that this model will be run. This example ranges from 1 to 24
ndays = 5
# this is the starting date in a datetime
date = datetime.datetime (2009, 11, 25, 0)
# this the time between outputs in the correct format (seconds)
tseas = 4*3600
# this is the time units. However I dont understand the date. Why it is 1970?
time_units = 'seconds since 1970-01-01'
# sets a smaller limit than between model outputs for when to force interpolation if it hasn't already occured
nsteps = 5
# controls the sampling frequency of the drifter tracks
N = 4
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

num_layers = 30
zpar = num_layers-1

# nwgom is a built i projection that is used. Also - galveston and nwgom-pyproj
proj = tracpy.tools.make_proj('nwgom')
# read in the grid
grid = tracpy.inout.readgrid (loc, proj, usespherical=True)
# initialization of the tracpy class
tp = Tracpy (loc, grid, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, time_units=time_units)
# input startinf locations in the long and latitude
lon0, lat0 = np.meshgrid(np.linspace(-98.5,-87.5,5), np.linspace(22.5,31,4))
# removes points not found within the domain
lon0, lat0 = tracpy.tools.check_points(lon0, lat0, tp.grid)
# timing that the grid was already read in
lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)
# plot the tracks
tracpy.plotting.tracks(lonp, latp, tp.name, tp.grid)

tracpy.plotting.hist(lonp, latp, tp.name, grid=tp.grid, which='hexbin', bins=(50,50))