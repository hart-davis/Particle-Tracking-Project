# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 09:01:46 2016

@author: support
"""
import os
import netCDF4 
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.cm
import matplotlib.pyplot as plt
from matplotlib import colors as c
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap
import cmocean as cm

nc = netCDF4.Dataset('/home/support/Downloads/20140101060000-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v02.0-fv01.0.nc')

x = nc.variables['lon'][:]
y = nc.variables['lat'][:]

m = Basemap(projection='cea',  # 2
                  urcrnrlon=35, 
                  llcrnrlon=17,
                  urcrnrlat=-23.2, 
                  llcrnrlat=-41.2, 
                  resolution='h')

x_map, y_map = m(*np.meshgrid(x, y))  # 3

#m.pcolormesh(x_map, y_map,  cmap = 'rainbow', vmin=-0.1, vmax=0.5)

m.drawcoastlines(linewidth=1.0)
m.fillcontinents()
m.drawcountries()
#m.drawparallels(np.arange(-40, -23, 1))
#m.drawmeridians(np.arange(17., 35, 1))
plt.colorbar()
plt.show()