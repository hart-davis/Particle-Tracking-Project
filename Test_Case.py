# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:21:38 2016

@author: support
"""
import os
import netCDF4 
import numpy as np
from numpy import where, nan
import pandas as pd
import xray
import matplotlib.cm
import matplotlib.pyplot as plt
from matplotlib import colors as c
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap
import cmocean as cm
import glob
from matplotlib import animation

#nc = netCDF4.Dataset('/home/support/Downloads/20140101060000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc')
#ncx = xray.open_dataset('/home/support/Downloads/20140101060000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc', decode_cf=False)

files = sorted(glob.glob('/home/support/Desktop/GlobCurrent_Data/20080101??????-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v02.0-fv01.0.nc'))
ds = xray.concat([xray.open_dataset(f) for f in files], 'time')


lat = ds.lat
lon = ds.lon
time = ds.time
u = ds.eastward_eulerian_current_velocity
v = ds.northward_eulerian_current_velocity

#lat = ncx.variables['lat'][:]
#lon = ncx.variables['lon'][:]

#time = ncx.variables['eastward_ekman_current_velocity']


#fig, ax = plt.subplots(figsize=[12,8])
#plt.pcolormesh(lon, lat, time, cmap='ocean', vmin=-1, vmax=1)

#plt.colorbar()
#plt.show()

#x = ncx.variables['lon'][:]
#y = ncx.variables['lat'][:]

#u= ncx.variables['eastward_ekman_current_velocity'][0]
#v = ncx.variables['northward_ekman_current_velocity'][0]

#u[where (u==ncx.eastward_ekman_current_velocity._FillValue)] = np.nan
#v[where (v==ncx.northward_ekman_current_velocity._FillValue)] = np.nan
#
velocity = np.sqrt(u**2+v**2)
#
m = Basemap(projection='cea',  # 2
                  urcrnrlon=35, 
                  llcrnrlon=17,
                  urcrnrlat=-23.2, 
                  llcrnrlat=-41.2, 
                  resolution='h')

x_map, y_map = m(*np.meshgrid(lon, lat))  # 3

for i in range(0,len(time)):
    m.pcolormesh(x_map, y_map, velocity[i], cmap = 'rainbow', vmin=0., vmax=2)
    m.drawcoastlines(linewidth=1.0)
    m.fillcontinents()

#ani = animation.FuncAnimation(x_map, y_map, velocity, interval=1000, blit=False)

#m.drawcountries()
#m.drawparallels(np.arange(-40, -23, 1))
#m.drawmeridians(np.arange(17., 35, 1))
    
    
    plt.colorbar()
    plt.savefig("day1_{:0}.jpeg".format(i))
    plt.clf()
    plt.show()
    

    
    

