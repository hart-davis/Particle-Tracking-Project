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

matplotlib.rcParams.update({'font.size': 10})
# this is the location of a high resolution model TXLA (texas something) Change in all files
loc = '/home/hart-davis/Desktop/Globcurrent_data_final.nc'
# number of days that this model will be run. This example ranges from 1 to 24
ndays = 1 
# this is the starting date in a datetime
date = datetime.datetime (2008, 1, 1, 0)
# this the time between outputs in the correct format (seconds)
tseas = 4*3600

time_units = 'seconds since 1970-01-01'
# sets a smaller limit than between model outputs for when to force interpolation if it hasn't already occured
nsteps = 5
# controls the sampling frequency of the drifter tracks
N = 4 #4
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

files = sorted(glob.glob('/home/hart-davis/Desktop/Glob_Current/2008010???????-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v02.0-fv01.0.nc'))

ds = xray.open_dataset(files[0]) 
lat = ds.lat.data
lon = ds.lon.data

west = 17 #360 #check the long and lati co-ordinates
east = 35 #378
south = -41.2
north = -23.2

ilon = where((lon >= west) & (lon <= east))
ilat = where((lat >= south) & (lat <= north))
ilon = squeeze(array(ilon))
ilat = squeeze(array(ilat))
lat = lat[ilat]
lon = lon[ilon]



# nwgom is a built i projection that is used. Also - galveston and nwgom-pyproj
proj = tracpy.tools.make_proj('cape_agulhas')
# read in the grid
grid = tracpy.inout.readgrid(loc, proj, usespherical = False)

x = lon
y = lat

xu = op.resize(x,0)
yu = y
xv = x
yv = op.resize(y,0)
XU, YU = np.meshgrid(xu, yu)
XV, YV = np.meshgrid(xv, yv)
XPSI, YPSI = np.meshgrid(xv, yu)
X, Y = np.meshgrid(x, y)
   
grid.x_u = np.asfortranarray(XU.T); grid.y_u = np.asfortranarray(YU.T)
grid.x_v = np.asfortranarray(XV.T); grid.y_v = np.asfortranarray(YV.T)
grid.x_rho = np.asfortranarray(X.T); grid.y_rho = np.asfortranarray(Y.T)
grid.x_psi = np.asfortranarray(XPSI.T); grid.y_psi = np.asfortranarray(YPSI.T)

grid.lon_rho, grid.lat_rho = grid.x_rho, grid.y_rho
grid.lon_u, grid.lat_u = grid.x_u, grid.y_u
grid.lon_v, grid.lat_v = grid.x_v, grid.y_v
grid.lon_psi, grid.lat_psi = grid.x_psi, grid.y_psi
grid.lon_rho = np.asfortranarray(grid.lon_rho)
grid.lat_rho = np.asfortranarray(grid.lat_rho)
grid.lon_u = np.asfortranarray(grid.lon_u)
grid.lat_u = np.asfortranarray(grid.lat_u)
grid.lon_v = np.asfortranarray(grid.lon_v)
grid.lat_v = np.asfortranarray(grid.lat_v)
grid.lon_psi = np.asfortranarray(grid.lon_psi)
grid.lat_psi = np.asfortranarray(grid.lat_psi)

grid.Cs_r = np.array([-1., 0.])
grid.sc_r = np.array([-1., 0.]) 

grid.h = np.asfortranarray(1.*np.ones(grid.x_rho.shape))
grid.zrt0 = np.asfortranarray(-0.5*np.ones((grid.x_rho.shape[0], grid.x_rho.shape[1], grid.Cs_r.size-1)))
grid.zwt0 = np.asfortranarray(-0.5*np.ones((grid.x_rho.shape[0], grid.x_rho.shape[1], grid.Cs_r.size)))
grid.dzt0 = np.asfortranarray(1.*np.ones((grid.x_rho.shape[0], grid.x_rho.shape[1], grid.Cs_r.size-1)))
 
grid.dxv = np.asfortranarray((x[1]-x[0])*np.ones(grid.x_v.shape))
grid.dyu = np.asfortranarray((y[1]-y[0])*np.ones(grid.y_u.shape))
grid.dxdy = np.asfortranarray((x[1]-x[0])*(y[1]-y[0])*np.ones(grid.y_rho.shape))
grid.imt = grid.x_rho.shape[0]
grid.jmt = grid.x_rho.shape[1]
grid.km = 1
grid.kmt = np.ones((grid.imt,grid.jmt),order='f')*grid.km

grid.X, grid.Y = x, y
grid.X = np.asfortranarray(grid.X)
grid.Y = np.asfortranarray(grid.Y)

grid.Vtransform = 1
grid.Vstretching = 1
grid.theta_s = 1
grid.theta_b = 0
grid.hc = 0


# initialization of the tracpy class
tp = Tracpy (loc, grid, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=0, doturb=0, do3d=0, z0=z0, zpar=0, time_units=time_units)
# input starting locations in the long and latitude
lon0, lat0 = np.meshgrid(np.linspace(15,35,20), np.linspace(-25,-35,15))
# removes points not found within the domain
lon0, lat0 = tracpy.tools.check_points(lon0, lat0, tp.grid)
# timing that the grid was already read in
lonp, latp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)
# plot the tracks
tracpy.plotting.tracks(lonp, latp, tp.name, tp.grid)

tracpy.plotting.hist(lonp, latp, tp.name, grid=tp.grid, which='hexbin', bins=(50,50))
