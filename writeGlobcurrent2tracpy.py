# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 12:13:25 2016

@author: hart-davis
"""

import xray
from numpy import where, squeeze, array, sqrt, zeros
import os
import glob
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from numpy import empty, size, unique, diff, asfortranarray
from datetime import datetime, timedelta
from tracpy import op
#files = sorted(glob.glob('netcdf/2008??????????-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v02.0-fv01.0.nc'))
files = sorted(glob.glob('/home/hart-davis/Desktop/Glob_Current/200801????????-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v02.0-fv01.0.nc'))

ds = xray.open_dataset(files[0]) # this might need tweeking 
lat = ds.lat.data
lon = ds.lon.data

west = -50 #360 #check the long and lati co-ordinates
east = 27 #378
south = -34
north = -33.5

ilon = where((lon >= west) & (lon <= east))
ilat = where((lat >= south) & (lat <= north))
ilon = squeeze(array(ilon))
ilat = squeeze(array(ilat))
lat = lat[ilat]
lon = lon[ilon]
ilon_u = ilon[:-1]
ilat_u = ilat[:]
ilon_v = ilon[:]
ilat_v = ilat[:-1]


# specify some variables to use to create dimensions in netcdf file
# eta = lat dimension
# xi = lon dimension
dim_eta_rho = len(lat)
dim_xi_rho = len(lon)
dim_eta_psi = len(lat)-1
dim_xi_psi = len(lon)-1
dim_eta_u = len(lat)
dim_xi_u = len(lon)-1
dim_eta_v = len(lat)-1
dim_xi_v = len(lon)
dim_ocean_time = len(files)
dim_s_rho = 1

# define the coriolis force
from numpy import pi
f = 4 * pi * np.sin(pi*lat/180) / (24*3600)
dims = [dim_eta_rho, dim_xi_rho]
cf = zeros(dims)
for i in range(0,dim_eta_rho):
	cf[:,i] = f

x = lon # at the rho point
y = lat # at the rho point
xu = op.resize(x,0)
yu = y
xv = x
yv = op.resize(y,0)
XU, YU = np.meshgrid(xu, yu)
XV, YV = np.meshgrid(xv, yv)
XPSI, YPSI = np.meshgrid(xv, yu)
X, Y = np.meshgrid(x, y)
   
x_u1 = np.asfortranarray(XU); y_u1 = np.asfortranarray(YU)
x_v1 = np.asfortranarray(XV); y_v1 = np.asfortranarray(YV)
x_rho1 = np.asfortranarray(X); y_rho1 = np.asfortranarray(Y)
x_psi1 = np.asfortranarray(XPSI); y_psi1 = np.asfortranarray(YPSI)

lon_rho1, lat_rho1 = x_rho1, y_rho1
lon_u1, lat_u1 = x_u1, y_u1
lon_v1, lat_v1 = x_v1, y_v1
lon_psi1, lat_psi1 = x_psi1, y_psi1
lon_rho1 = np.asfortranarray(lon_rho1)
lat_rho1 = np.asfortranarray(lat_rho1)
lon_u1 = np.asfortranarray(lon_u1)
lat_u1 = np.asfortranarray(lat_u1)
lon_v1 = np.asfortranarray(lon_v1)
lat_v1 = np.asfortranarray(lat_v1)
lon_psi1 = np.asfortranarray(lon_psi1)
lat_psi1 = np.asfortranarray(lat_psi1)

Cs_r1 = np.array([-1., 0.])
sc_r1 = np.array([-1., 0.]) 

h1 = np.asfortranarray(1.*np.ones(x_rho1.shape))
zrt01 = np.asfortranarray(-0.5*np.ones((x_rho1.shape[0], x_rho1.shape[1], Cs_r1.size-1)))
zwt01 = np.asfortranarray(-0.5*np.ones((x_rho1.shape[0], x_rho1.shape[1], Cs_r1.size)))
dzt01 = np.asfortranarray(1.*np.ones((x_rho1.shape[0], x_rho1.shape[1], Cs_r1.size-1)))
 
dxv1 = np.asfortranarray((x[1]-x[0])*np.ones(x_v1.shape))
dyu1 = np.asfortranarray((y[1]-y[0])*np.ones(y_u1.shape))
dxdy1 = np.asfortranarray((x[1]-x[0])*(y[1]-y[0])*np.ones(y_rho1.shape))
imt1 = x_rho1.shape[0]
jmt1 = x_rho1.shape[1]
km1 = 1
kmt1 = np.ones((imt1,jmt1),order='f')*km1

X1, Y1 = np.meshgrid(x, y)
X1 = np.asfortranarray(X1)
Y1 = np.asfortranarray(Y1)

Vtransform1 = 1
Vstretching1 = 1
theta_s1 = 1
theta_b1 = 0
hc1 = 0

dimsPM = [dim_eta_rho, dim_xi_rho]
pm1 = zeros(dimsPM)
pm1[:] = 0.25

dimsPN = [dim_eta_rho, dim_xi_rho]
pn1 = zeros(dimsPN)
pn1[:] = 0.25

time = zeros(dim_ocean_time)
u1 = zeros([dim_ocean_time, dim_eta_u, dim_xi_u])
v1 = zeros([dim_ocean_time, dim_eta_v, dim_xi_v])
for f in range(0,len(files)):
	ds = xray.open_dataset(files[f])
	#print files[f]
	time[f] = (ds.time.data - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
	u1[f,:,:] = ds.eastward_eulerian_current_velocity[:,ilat_u,ilon_u]
	v1[f,:,:] = ds.northward_eulerian_current_velocity[:,ilat_v,ilon_v]


z1 = zeros(dims)
z1=z1.astype(int) 

levels = [0]

fname = 'Port_Elizabeth.nc'
rootgrp = Dataset(fname, "w", format="NETCDF4")

# define dimensions
eta_rho = rootgrp.createDimension("eta_rho", dim_eta_rho)
xi_rho = rootgrp.createDimension("xi_rho", dim_xi_rho)
eta_psi = rootgrp.createDimension("eta_psi", dim_eta_psi)
xi_psi = rootgrp.createDimension("xi_psi", dim_xi_psi)
eta_u = rootgrp.createDimension("eta_u", dim_eta_u)
xi_u = rootgrp.createDimension("xi_u", dim_xi_u)
eta_v = rootgrp.createDimension("eta_v", dim_eta_v)
xi_v = rootgrp.createDimension("xi_v", dim_xi_v)
ocean_time = rootgrp.createDimension("ocean_time", dim_ocean_time)
s_rho = rootgrp.createDimension("s_rho", dim_s_rho)

# create variables
lon_psi = rootgrp.createVariable("lon_psi","f8",("eta_psi","xi_psi",))
lat_psi = rootgrp.createVariable("lat_psi","f8",("eta_psi","xi_psi",))
x_psi = rootgrp.createVariable("x_psi","f8",("eta_psi","xi_psi",))
y_psi = rootgrp.createVariable("y_psi","f8",("eta_psi","xi_psi",))
sc_r = rootgrp.createVariable("sc_r","f8",("s_rho",))
Cs_r = rootgrp.createVariable("Cs_r","f8",("s_rho",))
zrt0 = rootgrp.createVariable("zrt0","f8",("eta_rho", "xi_rho", "s_rho",))
zwt0 = rootgrp.createVariable("zwt0","f8",("eta_rho", "xi_rho", "s_rho",))
dzt0 = rootgrp.createVariable("dzt0","f8",("eta_rho", "xi_rho", "s_rho",))
dxv = rootgrp.createVariable("dxv","f8",("eta_v", "xi_v",))
dyu = rootgrp.createVariable("dyu","f8",("eta_u","xi_u",))
dxdy = rootgrp.createVariable("dxdy","f8",("eta_rho","xi_rho"))
kmt = rootgrp.createVariable("kmt","f8",("eta_rho","xi_rho",))
km = rootgrp.createVariable("km","f8",)
jmt = rootgrp.createVariable("jmt","f8",("eta_rho",))
imt = rootgrp.createVariable("imt","f8",("xi_rho",))
Vtransform = rootgrp.createVariable("Vtransform","f8",("s_rho",))
Vstretching = rootgrp.createVariable("Vstretching","f8",("s_rho",))
theta_s = rootgrp.createVariable("theta_s","f8",("s_rho"))
theta_b = rootgrp.createVariable("theta_b","f8",("s_rho"))
lon_v = rootgrp.createVariable("lon_v","f8",("eta_v", "xi_v",))
lat_v = rootgrp.createVariable("lat_v","f8",("eta_v", "xi_v",))
lon_u = rootgrp.createVariable("lon_u","f8",("eta_u", "xi_u",))
lat_u = rootgrp.createVariable("lat_u","f8",("eta_u", "xi_u",))
f = rootgrp.createVariable("f","f8",("eta_rho","xi_rho",))
ocean_time = rootgrp.createVariable("ocean_time","f8",("ocean_time",))
u = rootgrp.createVariable("u","f8",("ocean_time","s_rho","eta_u","xi_u",),fill_value = 0)
v = rootgrp.createVariable("v","f8",("ocean_time","s_rho","eta_v","xi_v",),fill_value = 0)
h = rootgrp.createVariable("h", "f8", ("eta_rho","xi_rho",))
hc = rootgrp.createVariable("hc", "f8")
s_rho = rootgrp.createVariable("s_rho","f8",("s_rho",))
lat_rho = rootgrp.createVariable("lat_rho","f8",("eta_rho","xi_rho",))
lon_rho = rootgrp.createVariable("lon_rho","f8",("eta_rho","xi_rho",))
pm = rootgrp.createVariable("pm","f8",("eta_rho","xi_rho",))
#mask_psi = rootgrp.createVariable("mask_psi", "f8", ("mask_psi",))
#mask_rho = rootgrp.createVariable("mask_rho", "f8", ("mask_rho",))
#mask_u = rootgrp.createVariable("mask_u", "f8", ("mask_u",))
#mask_v = rootgrp.createVariable("mask_v", "f8", ("mask_v",))
pn = rootgrp.createVariable("pn","f8",("eta_rho","xi_rho",))
x_rho = rootgrp.createVariable("x_rho","f8",("eta_rho","xi_rho",))
y_rho = rootgrp.createVariable("y_rho","f8",("eta_rho","xi_rho",))
x_u = rootgrp.createVariable("x_u","f8",("eta_u", "xi_u",))
y_u = rootgrp.createVariable("y_u","f8",("eta_u", "xi_u",))
x_v = rootgrp.createVariable("x_v","f8",("eta_v", "xi_v",))
y_v = rootgrp.createVariable("y_v","f8",("eta_v", "xi_v",))
#angle = rootgrp.createVariable("angle","f8",("lat_rho","lon_rho",))
##lon_vert = rootgrp.createVariable("x_vert","f8",("lat","lon"))
##lat_vert = rootgrp.createVariable("y_vert","f8",("lat","lon"))
zp = rootgrp.createVariable("zp", "f8",("eta_rho","xi_rho"))
z = rootgrp.createVariable("z", "f8",("ocean_time","s_rho","eta_rho","xi_rho",),fill_value = 9.999999933815813e+36)

import time as t
# global
rootgrp.description = "Mikes Globcurrent subset + Grid Data"
rootgrp.history = "Created " + t.ctime(t.time())
rootgrp.time_coverage_start = "1-Jan-2008"
rootgrp.time_coverage_end = "31-Aug-2008"
rootgrp.spatialResolution = str(float(ds.spatial_resolution)) #RESOLUTION OF GLOBCURRENT!

# by variable
# longitudes
lon_v.units = "degrees_east"
lon_v.long_name = "longitude of V-points"
lon_v.standard_name = "longitude"
lon_v.field = "lon_v, scalar"
lon_u.units = "degrees_east"
lon_u.long_name = "longitude of U-points"
lon_u.standard_name = "longitude"
lon_u.field = "lon_u, scalar"
# latitudes
lat_v.units = "degrees_north"
lat_v.long_name = "latitude of V-points"
lat_v.standard_name = "latitude"
lat_v.field = "lat_v, scalar"
lat_u.units = "degrees_north"
lat_u.long_name = "latitude of U-points"
lat_u.standard_name = "latitude"
lat_u.field = "lat_u, scalar"

lat_rho.units = "degrees_north"
lat_rho.long_name = "latitude of RHO-points"
lat_rho.standard_name = "latitude"
lat_rho.field = "lat_rho, scalar"

lon_rho.units = "degrees_east"
lon_rho.long_name = "longitude of RHO-points"
lon_rho.standard_name = "longitude"
lon_rho.field = "lon_rho, scalar"
# times
ocean_time.long_name = "time since initialization"
ocean_time.units = "seconds since 1970-01-01 00:00:00"
ocean_time.calendar = "gregorian"
ocean_time.field = "time, scalar, series"
ocean_time.standard_name = "time"

# sst
u.long_name = "u-momentum component"
u.units = "meter second-1"
u.time = "ocean_time"
u.coordinates = "ocean_time s_rho lat_u lon_u"
u.field = "u-velocity, scalar, series"
u.fill_value = 9.999999933815813e+36

v.long_name = "v-momentum component"
v.units = "meter second-1"
v.time = "ocean_time"
v.coordinates = "ocean_time s_rho lat_v lon_v"
v.field = "v-velocity, scalar, series"
v.fill_value = 9.999999933815813e+36

s_rho.long_name = "S-coordinate at RHO-points"
s_rho.valid_min = "-1.0"
s_rho.valid_max = "0.0"
s_rho.positive = "up"
s_rho.standard_name = "ocean_s_coordinate_g1"
s_rho.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"
s_rho.field = "s_rho, scalar"

pm.long_name = "curvilinear coordinate metric in XI"
pm.units = "met,er-1"
pm.coordinates = "lat_rho lon_rho"
pm.field = "pm, scalar"

f.long_name = "Coriolis parameter at RHO-points"
f.units = "second-1"
f.coordinates = "lat_rho lon_rho"
f.field = "coriolis,scalar"

h.long_name = "bathymetry at RHO-points"
h.units = "meter"
h.coordinates = "lat_rho lon_rho"
h.field = "bath,scalar"

#mask_psi.long_name = "mask on psi-points"
#mask_psi.flag_value = "0.0, 1.0"
#mask_psi.flag_meanings = "land water"
#mask_psi.coordinates = "lat_psi lon_psi"
#
#mask_rho.long_name = "mask on rho-points"
#mask_rho.flag_value = "0.0, 1.0"
#mask_rho.flag_meanings = "land water"
#mask_rho.coordinates = "lat_rho lon_rho"
#
#mask_u.long_name = "mask on u-points"
#mask_u.flag_value = "0.0, 1.0"
#mask_u.flag_meanings = "land water"
#mask_u.coordinates = "lat_u lon_u"
#
#mask_v.long_name = "mask on v-points"
#mask_v.flag_value = "0.0, 1.0"
#mask_v.flag_meanings = "land water"
#mask_v.coordinates = "lat_v lon_v"

pn.long_name = "curvilinear coordinate metric in ETA"
pn.units = "meter-1"
pn.coordinates = "lat_rho lon_rho"
pn.field = "pn, scalar"

x_rho.long_name = "x_rho"
x_rho.units = "degrees"
x_rho.coordinates = "x_rho lat_rho lon_rho"
x_rho.field = "x_rho lat_rho lon_rho"

y_rho.long_name = "y_rho"
y_rho.units = "degrees"
y_rho.coordinates = "y_rho lat lon"
y_rho.field = "y_rho lat lon"

#xr.long_name = "xr"
#xr.units = "degrees"
#xr.coordinates = "xr lat lon"
#xr.field = "xr lat lon"
#
#yr.long_name = "yr"
#yr.units = "degrees"
#yr.coordinates = "yr lat lon"
#yr.field = "yr lat lon"

x_u.long_name = "x_u"
x_u.units = "degrees"
x_u.coordinates = "x_u lat_u lon_u"
x_u.field = "x_u lat_u lon_u"

y_u.long_name = "y_u"
y_u.units = "degrees"
y_u.coordinates = "y_u lat_u lon_u"
y_u.field = "y_u lat_u lon_u"

x_v.long_name = "x_v"
x_v.units = "degrees"
x_v.coordinates = "x_v lat_v lon_v"
x_v.field = "x_v lat_v lon_v"

y_v.long_name = "y_v"
y_v.units = "degrees"
y_v.coordinates = "y_v lat_v lon_v"
y_v.field = "y_v lat_v lon_v"

#angle.long_name = "angle between XI-axis and EAST"
#angle.units = "radians"
#angle.coordinates = "lat_rho lon_rho"
#angle.field = "angle, scalar"

#lon_vert.long_name = "x vert"
#lon_vert.units = "0"
#lon_vert.coordinates = "lat lon"
#lon_vert.field = "lat lon"

#lat_vert.long_name = "y vert"
#lat_vert.units = "0"
#lat_vert.coordinates = "lat lon"
#lat_vert.field = "lat lon"

z.long_name = "vertical velocity"
z.units = "m.s-1"
z.coordinates = "ocean_time s_rho lat_rho lon_rho"
v.field = "v-velocity, scalar, series"
v.fill_value = 9.999999933815813e+36

Vtransform.long_name = "vertical terrain-following transformation equation"

Vstretching.long_name = "vertical terrain-following stretching function"

# write data to variables
#lon[0:len(lon)] = lonocean_time[0:len(time)] = time
ocean_time[0:len(time)] = time
u[0:dim_ocean_time, 0:dim_s_rho, 0:dim_eta_u, 0:dim_xi_u] = u1
v[0:dim_ocean_time, 0:dim_s_rho, 0:dim_eta_v, 0:dim_xi_v] = v1
s_rho[0:dim_s_rho] = levels
pm[0:dim_eta_rho, 0:dim_xi_rho] = pm1
f[0:dim_eta_rho, 0:dim_xi_rho] = cf
km = km1
lon_psi[0:dim_eta_psi, 0:dim_xi_psi] = lon_psi1[0:dim_eta_psi, 0:dim_xi_psi]
lat_psi[0:dim_eta_psi, 0:dim_xi_psi] = lat_psi1[0:dim_eta_psi, 0:dim_xi_psi]
x_psi[0:dim_eta_psi, 0:dim_xi_psi] = x_psi1[0:dim_eta_psi, 0:dim_xi_psi]
y_psi[0:dim_eta_psi, 0:dim_xi_psi] = y_psi1[0:dim_eta_psi, 0:dim_xi_psi]
sc_r[0:dim_s_rho] = sc_r1[0:dim_s_rho]
Cs_r[0:dim_s_rho] = Cs_r1[0:dim_s_rho]
# dodgy as fuck
zrt0[0:dim_eta_rho, 0:dim_xi_rho, 0] = zrt01
zwt0 = zwt01
dzt0[0:dim_eta_rho, 0:dim_xi_rho, 0] = dzt01
dxv[0:dim_eta_psi, 0:dim_xi_psi] = dxv1[0:dim_eta_psi, 0:dim_xi_psi]
dxdy[:, :] = dxdy1[:, :]
kmt = kmt1[:,:]
jmt = jmt1
imt = imt1
Vtransform = Vtransform1
Vstretching= Vstretching1
theta_s = theta_s1
theta_b = theta_b1
lon_v[0:dim_eta_v, 0:dim_xi_v] = lon_v1
lat_v[0:dim_eta_v, 0:dim_xi_v] = lat_v1
lon_u[0:dim_eta_u, 0:dim_xi_u] = lon_u1
lat_u[0:dim_eta_u, 0:dim_xi_u] = lat_u1

h[:, :] = h1[:, :]
hc = hc1
s_rho = levels
lat_rho[0:dim_eta_rho, 0:dim_xi_rho] = lat_rho1[0:dim_eta_rho, 0:dim_xi_rho]
lon_rho[0:dim_eta_rho, 0:dim_xi_rho] = lon_rho1[0:dim_eta_rho, 0:dim_xi_rho]
pn[0:dim_eta_rho, 0:dim_xi_rho] = pn1[0:dim_eta_rho, 0:dim_xi_rho]
x_rho[0:dim_eta_rho, 0:dim_xi_rho] = x_rho1[0:dim_eta_rho, 0:dim_xi_rho]
y_rho[0:dim_eta_rho, 0:dim_xi_rho] = y_rho1[0:dim_eta_rho, 0:dim_xi_rho]
x_u[0:dim_eta_u, 0:dim_xi_u] = x_u1
y_u[0:dim_eta_u, 0:dim_xi_u] = y_u1
x_v[0:dim_eta_v, 0:dim_xi_v] = x_v1
y_v[0:dim_eta_v, 0:dim_xi_v] = y_v1
dyu= dyu1[:, :]

z = z1
#zp[zp1] = zp1

rootgrp.close()
