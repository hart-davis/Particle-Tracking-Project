import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature

def drifter_plotting3d(x,y,c,ctype=None,style=style,
                       cmap=plt.cm.Spectral_r,level=None,vmin=None,vmax=None,
                       bins=None,res='110m',proj=ccrs.Mercator(),cbar=False):
    
    """
    x : longitudinal coordinate of  point/s
    y : latitudinal  coordinate of point/s
    c : property of coordinate. e.g. velocity at x,y.

    ctype: is the structure of x,y and c. If x,y and c are one-d arrays then ctype = None. [default] 
              If x,y,c are in 2-d arrays, then ctype = two_dim
    
    style: default = None, which is a scatter plot. Options are scatter, pcolormesh, contourf. This can be 
           developed and will be upgraded. 
    
    cmap: color map for c. Default is plt.cm.Spectral_r
    
    level: For contour plots for the number of contours to be used. 
    
    vmin/vmax: Minimum and maximum values to illustrate in plots. 
    
    bins: The bin value for pcolormesh.
    
    res: Default is 110m. Can be set as 10m ,50n and 110m.
    
    proj: Projection to be used in cartopy. Options available in cartopy. 
          Land and States are plotted by default.
              """
    print "Plot Set as : "
    land = cfeature.NaturalEarthFeature('physical', 'land', res,
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    states = cfeature.NaturalEarthFeature('cultural',
                                              name='admin_1_states_provinces_lines',
                                              scale=res, facecolor='none')
    proj=proj
    res =res

    fig = plt.figure(figsize=(20, 7.7))
    ax = plt.subplot(111, projection=proj)
    pc = ccrs.PlateCarree()

    extent=[0, 48, -48, -20]
    ax.set_extent(extent, pc)

    ax.add_feature(land, facecolor='slategray',zorder=4)
    ax.coastlines(resolution=res)  
    c = c[:len(x)]
    if style == None:
        print "None"
        ax.plot(   x,y, marker='.', ms=0.2, c='k', lw=0,transform=pc)
        ax.plot(   x,y, marker='.', ms=0.1, c='blue', lw=0,transform=pc)
        
    elif style == scatter:
        ax.scatter(x,y,c=c,s=1,transform=pc,cmap=cmap,vmin=vmin,vmax=vmax)
    elif style == 'pcolormesh':  
        print "Pcolormesh"
        if ctype == None:
            if  bins == None:
                gridxt = np.linspace(x.min(), x.max(), (x.max()-x.min())/2)
                gridyt = np.linspace(y.min(), y.max(), (y.max()-y.min())/2)
            elif bins== bins:
                gridxt = np.linspace(x.min(), x.max(), bins[0])
                gridyt = np.linspace(y.min(), y.max(), bins[1])
                
            longitudes, latitudes = np.meshgrid(gridxt, gridyt, sparse=True)
            #longitudes = np.transpose(longitudes)
            longitudes = np.squeeze(longitudes)
            latitudes  =  np.squeeze(latitudes)
            indxcheck  =  np.zeros([len(longitudes),len(latitudes)])
            vel   =   np.zeros([len(longitudes),len(latitudes)])
            # loop through the grid, find lon/lat values for each 0.25 grid cell and average corresponding velocities
            for i in     range(0, len(longitudes)-1, 1):
                for j in range(0, len(latitudes)-1, 1):
                    indx = np.where(  (x[:len(c)]>=longitudes[i]) 
                                    & (x[:len(c)]<longitudes[i+1]) 
                                    & (y[:len(c)]>=latitudes[j]) 
                                    & (y[:len(c)]<latitudes[j+1]))[0]

                    indxcheck[i,j]=len(indx)
                    # if no drifters found put a nan in the grid cell
                    if indx.size == 0:
                        vel[i,j] = 0
                    else:
                        vel[i,j] = np.mean(velocity[indx])
            velo = np.transpose(vel)
            velo = ma.masked_where(np.isnan(velo),velo)
            ss= ax.pcolormesh(gridxt,gridyt,velo,transform=pc,cmap=cmap,vmin=vmin,vmax=vmax)
        elif ctype == two_dim:
            ss= ax.pcolormesh(x,y,c,transform=pc,cmap=cmap,vmin=vmin,vmax=vmax)
    elif style == 'contour':
        print "Contour"
        if ctype == None:
            
            gridxt = np.linspace(x.min(), x.max(), (x.max()-x.min())/2)
            gridyt = np.linspace(y.min(), y.max(), (y.max()-y.min())/2)
            longitudes, latitudes = np.meshgrid(gridxt, gridyt, sparse=True)
            #longitudes = np.transpose(longitudes)
            longitudes = np.squeeze(longitudes)
            latitudes  =  np.squeeze(latitudes)
            indxcheck  =  np.zeros([len(longitudes),len(latitudes)])
            vel   =   np.zeros([len(longitudes),len(latitudes)])
            # loop through the grid, find lon/lat values for each 0.25 grid cell and average corresponding velocities
            for i in     range(0, len(longitudes)-1, 1):
                for j in range(0, len(latitudes)-1, 1):
                    indx = np.where(  (x[:len(c)]>=longitudes[i]) 
                                    & (x[:len(c)]<longitudes[i+1]) 
                                    & (y[:len(c)]>=latitudes[j]) 
                                    & (y[:len(c)]<latitudes[j+1]))[0]

                    indxcheck[i,j]=len(indx)
                    # if no drifters found put a nan in the grid cell
                    if indx.size == 0:
                        vel[i,j] = 0
                    else:
                        vel[i,j] = np.mean(velocity[indx])
            velo = np.sqrt(vel**2)
            velo = np.transpose(velo)
            velo = ma.masked_where(np.isnan(velo),velo)
            ss= ax.contourf(gridxt,gridyt,velo,transform=pc,cmap=cmap,vmin=vmin,vmax=vmax,
                            levels=np.linspace(vmin,vmax,level))
            ass = ax.contour(gridxt,gridyt,velo,transform=pc,colors='black',vmin=vmin,vmax=vmax,
                            levels=np.linspace(vmin,vmax,level))
        elif ctype == two_dim:
            ss= ax.contourf(x,y, c,transform=pc,cmap=cmap)
            ass = ax.contour(x,y,c,transform=pc,colors='black')
            
        if cbar == True:
            cbar=plt.colorbar(ss)
            cbar.vmax=vmax
