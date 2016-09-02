# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:21:38 2016

@author: support
"""
import os
import netCDF4 
import numpy as np
import pandas as pd
import xray
import matplotlib.cm
import matplotlib.pyplot as plt
from matplotlib import colors as c
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap
import cmocean as cm

nc = netCDF4.Dataset('/home/support/Downloads/20140101060000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc')
ncx = xray.open_dataset('/home/support/Downloads/20140101060000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc', decode_cf=False)
ncx

#lat = ncx.variables['lat'][:]
#lon = ncx.variables['lon'][:]

#time = ncx.variables['eastward_ekman_current_velocity']


#fig, ax = plt.subplots(figsize=[12,8])
#plt.pcolormesh(lon, lat, time, cmap='ocean', vmin=-1, vmax=1)

#plt.colorbar()
#plt.show()

x = ncx.variables['lon'][:]
y = ncx.variables['lat'][:]

east= ncx.variables['eastward_ekman_current_velocity'][0]

north = ncx.variables['northward_ekman_current_velocity'][0]

east_1 = np.sqrt(east**2)
north_1 = np.sqrt(north**2)

total = north_1 + east_1

m = Basemap(projection='cea',  # 2
                  urcrnrlon=35, 
                  llcrnrlon=17,
                  urcrnrlat=-23.2, 
                  llcrnrlat=-41.2, 
                  resolution='h')

x_map, y_map = m(*np.meshgrid(x, y))  # 3

m.pcolormesh(x_map, y_map, total, cmap = 'rainbow', vmin=0., vmax=0.1)

m.drawcoastlines(linewidth=1.0)
m.fillcontinents()
m.drawcountries()
#m.drawparallels(np.arange(-40, -23, 1))
#m.drawmeridians(np.arange(17., 35, 1))
plt.colorbar()
plt.show()

