from datetime import datetime
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import matplotlib.pyplot as plt
from matplotlib import patheffects
import metpy
from metpy.plots import colortables
import xarray as xr
from xarray.backends import NetCDF4DataStore
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
%matplotlib inline
# Product Name Reference Guide
# ACHA - Cloud Top Height: 'HT'    
# ACHT - Cloud Top Temperature: 'TEMP'
# ACTP - Cloud Top Phase: 'Phase'
# ADP - Aerosol Detection: 'Smoke'
# AOD - Aerosol Optical Depth: 'AOD'    
# COD - Cloud Optical Depth: 'COD'    
# CPS - Cloud Particle Size: 'PSD'
# CTP - Cloud Top Pressure: 'PRES'
# DMW - Derived Motion Winds: 'pressure','temperature', 'wind_direction', 'wind_speed'    
# DSI - Derived Stability Indices: 'CAPE', 'KI', 'LI', 'SI', 'TT' 
# FDC - Fire-Hot Spot Characterization: 'Area', 'Mask', 'Power', 'Temp'    
# FSC - Snow Cover: 'FSC'    
# LST - Land Surface (Skin) Temperature: 'LST'    
# RRQPE - Rainfall Rate - Quantitative Prediction Estimate: 'RRQPE'    
# RSR - Reflected Shortwave Radiation: 'RSR'    
# SST - Sea Surface (Skin) Temperature: 'SST'
# TPW - Total Precipitable Water: 'TPW'
# VAA - Volcanic Ash: 'VAH', 'VAML'    
def read_g16_product(product)
    domain='CONUS'
    # ACHAF - Cloud Top Height: 'HT'    
    elif product == "ACHA":
        longname = "CloudTopHeight"
        variable = 'HT'
        vmin = 0
        vmax = 15000
        cmap = "rainbow"

    # ACHTF - Cloud Top Temperature: 'TEMP'
    elif product == "ACHT":
        longname = "CloudTopTemperature"
        variable = 'TEMP' 
        vmin = 180
        vmax = 300
        cmap = "jet"
    
    # ACTPF - Cloud Top Phase: 'Phase'
    elif product == "ACTP":
        longname = "CloudTopPhase"
        variable = 'Phase' 
        vmin = 0
        vmax = 5
        cmap = "jet"     
    
    # ADPF - Aerosol Detection: 'Smoke'
    elif product == "ADP":
        longname = "AerosolDetection"
        variable = 'Smoke' 
        vmin = 0
        vmax = 255
        cmap = "jet"
    
    # AODF - Aerosol Optical Depth: 'AOD'    
    elif product == "AOD":
        longname = "AerosolOpticalDepth"
        variable = 'AOD'
        vmin = 0
        vmax = 2
        cmap = "rainbow"  
    
    # CODF - Cloud Optical Depth: 'COD'    
    elif product == "COD":
        longname = "CloudOpticalDepth"
        variable = 'COD'
        vmin = 0
        vmax = 100
        cmap = "jet"  
    
    # CPSF - Cloud Particle Size: 'PSD'
    elif product == "CPS":
        longname = "CloudParticleSize"
        variable = 'PSD'
        vmin = 0
        vmax = 80
        cmap = "rainbow"
    
    # CTPF - Cloud Top Pressure: 'PRES'
    elif product == "CTP":
        longname = "CloudTopPressure"
        variable = 'PRES'
        vmin = 0
        vmax = 1100
        cmap = "rainbow"
    
    # DMWF - Derived Motion Winds: 'pressure','temperature', 'wind_direction', 'wind_speed'    
    elif product == "DMW":
        longname = "DerivedWindMotion"
        variable = 'wind_direction' 
    
    # DSIF - Derived Stability Indices: 'CAPE', 'KI', 'LI', 'SI', 'TT' 
    elif product == "DSI":
        longname = "DerivedStabilityIndices"
        variable = 'CAPE' 
        vmin = 0
        vmax = 1000
        cmap = "jet" 
        
        #variable = 'KI'
        #vmin = -50
        #vmax = 50
        #cmap = "jet"
        
        #variable = 'LI'
        #vmin = -10
        #vmax = 30
        #cmap = "jet"
        
        #variable = 'SI'
        #vmin = -10
        #vmax = 25
        #cmap = "jet"
        
        #variable = 'TT'
        #vmin = -10
        #vmax = 60
        #cmap = "jet"
    
    # FDCF - Fire-Hot Spot Characterization: 'Area', 'Mask', 'Power', 'Temp'    
    elif product == "FDC":
        #variable = 'Area' 
        longname = "FireHotSpot"
        variable = 'Mask' 
        vmin = 0
        vmax = 255
        cmap = "jet"
        
        #variable = 'Power' 
        #variable = 'Temp' 
    
    # FSCF - Snow Cover: 'FSC'    
    elif product == "FSC":
        longname = "snow_cover_not_yet_available"
        variable = 'FSC' 
        vmin = 0
        vmax = 1
        cmap = "jet"
    
    # LSTF - Land Surface (Skin) Temperature: 'LST'    
    elif product == "LST":
        longname = "LandSurfaceTemperature"
        variable = 'LST'    
        vmin = 213
        vmax = 330
        cmap = "jet"
    
    # RRQPEF - Rainfall Rate - Quantitative Prediction Estimate: 'RRQPE'    
    elif product == "RRQPE":
        longname = "RainRateQPE"
        variable = 'RRQPE' 
        vmin = 0
        vmax = 35
        cmap = "jet"
    
    # RSR - Reflected Shortwave Radiation: 'RSR'    
    #elif product == "RSRF":
    #    variable = 'RSR'       
    
    # SSTF - Sea Surface (Skin) Temperature: 'SST'
    elif product == "SST":
        longname = "SeaSurfaceTemperature"
        variable = 'SST'
        vmin = 268
        vmax = 308
        cmap = "jet"
        domain='FullDisk'
    
    # TPWF - Total Precipitable Water: 'TPW'
    elif product == "TPW":
        longname = "TotalPrecipitableWater"
        variable = 'TPW' 
        vmin = 0
        vmax = 60
        cmap = "jet"
    
    # VAAF - Volcanic Ash: 'VAH', 'VAML'    
    elif product == "VAA":
        #variable = 'VAH'
        #vmin = 0
        #vmax = 20000
        #cmap = "jet" 
        longname = "VolcanicAshDetection"
        variable = 'VAML' 
        vmin = 0
        vmax = 100
        cmap = "jet"
        domain='FullDisk'
        
    nowdate = datetime.utcnow()
    cat = TDSCatalog('http://thredds-jumbo.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/Products/' + longname + '/' + domain + '/' + \
                      str(nowdate.year) + str("%02d"%nowdate.month) + str("%02d"%nowdate.day) + '/catalog.xml')
    dataset_name = sorted(cat.datasets.keys())[-1]
    dataset = cat.datasets[dataset_name]
    ds = dataset.remote_access(service='OPENDAP')
    ds = NetCDF4DataStore(ds)
    ds = xr.open_dataset(ds)
    
    
    dat = ds.metpy.parse_cf(variable)
    proj = dat.metpy.cartopy_crs
    
    dat = dat.where(dqf == 0)
    dat = dat.where(dat.variable > vmin)
    dat = dat.where(dat.variable < vmax)
    dat = dat - 273.15
    # Plot in Mercator
    newproj = ccrs.Mercator()
    
    if domain=='CONUS':
        toptext = 0.857
        toptextleft = 0.13
        toptextright = 0.7
        bottomtextleft = 0.13
        bottomtextheight = 0.155
        toprecx = 0.125
        toprecy = 0.755
        bottomrecx = 0.125
        bottomrecy = 0.25
        
        xoffmin = 40000
        xoffmax = -30000
        yoffmin = 55000
        yoffmax = -5000
    elif domain =='FullDisk':
        toptext = 0.857
        toptextleft = 0.13
        toptextright = 0.7
        bottomtextleft = 0.13
        bottomtextheight = 0.155
        toprecx = 0.125
        toprecy = 0.755
        bottomrecx = 0.125
        bottomrecy = 0.25
        xoffmin = 0
        xoffmax = 0
        yoffmin = 0
        yoffmax = 0
    
    if product=="SST":
        toptext = 0.857
        toptextleft = 0.13
        toptextright = 0.7
        bottomtextleft = 0.13
        bottomtextheight = 0.155
        toprecx = 0.125
        toprecy = 0.755
        bottomrecx = 0.125
        bottomrecy = 0.25
        
        xoffmin = 4000000
        xoffmax = -3200000
        yoffmin = 5500000
        yoffmax = -650000

    fig = plt.figure(figsize=[16, 12], dpi=100)
    ax = fig.add_subplot(1,1,1, projection=newproj)
    im = ax.pcolormesh(dat['x'], dat['y'], dat, cmap='jet', transform=proj, vmin=vmin, vmax=vmax)
    ax.set_extent((dat['x'].min() + xoffmin, dat['x'].max() + xoffmax, dat['y'].min() + yoffmin, dat['y'].max() + yoffmax), crs=proj)
    cbaxes = inset_axes(ax, width="100%", height="4%", loc='lower center', borderpad=0) 
    cb1 = fig.colorbar(im, orientation='horizontal', cax=cbaxes, ticks=[0, 5, 10, 15, 20, 25, 30, 35])
    cb1.ax.set_xticklabels(['0', '5', '10', '15', '20', '25', '30', '35'])
    cb1.outline.set_visible(False) # Remove the colorbar outline
    cb1.ax.tick_params(width = 0) # Remove the colorbar ticks 
    cb1.ax.xaxis.set_tick_params(pad=-15.5) # Put the colobar labels inside the colorbar
    ax.set_title("")
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, edgecolor='gray')
    ax.add_feature(cfeature.STATES, linestyle=':', edgecolor='gray')
    
    # top rectangle
    fig.patches.extend([plt.Rectangle((toprecx,toprecy),0.775,0.022,
                                  fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    # bottom rectangle
    fig.patches.extend([plt.Rectangle((bottomrecx,bottomrecy),0.775,0.022,
                                  fill=True, color='darkslateblue', alpha=1, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
    request = cimgt.GoogleTiles(url="https://cartodb-basemaps-d.global.ssl.fastly.net/dark_all/{z}/{x}/{y}.png")
    
    clabeltext='Radial Velocity [m/s]'
    title = 'CEMA Base Radial Velocity'
    timestr = "random time" #local.strftime('%Y-%m-%d %H:%M ') + et
    
    ax.add_image(request, 7, zorder=0, interpolation='none')
    fig.text(toptextleft, bottomtextheight,clabeltext,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    fig.text(toptextright, toptext,timestr,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    fig.text(bottomtextleft, toptext,title,horizontalalignment='left', color = 'white', size=10, zorder=2000)
    
    fig.savefig("Downloads/sst.png", dpi=100, bbox_inches='tight')







