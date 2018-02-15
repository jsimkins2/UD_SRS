%matplotlib inline
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
from pyproj import Proj     
from scipy import ndimage
from scipy import stats

# initial block mean allows us to downsample the Red Band which is a higher resolution than the others
def block_mean(ar, fact):
    assert isinstance(fact, int), type(fact)
    sx, sy = ar.shape
    X, Y = np.ogrid[0:sx, 0:sy]
    regions = sy/fact * (X/fact) + Y/fact
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res.shape = (sx/fact, sy/fact)
    return res





# C is for Conus File
C_file = 'Documents/Delaware/goesR/satdata/OR_ABI-L2-CMIPC-M3C02_G16_s20180431702207_e20180431704580_c20180431705095.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
R = C.variables['CMI'][:].data # Band 2 is red (0.64 um)
R = block_mean(R, 2)

C_file = 'Documents/Delaware/goesR/satdata/OR_ABI-L2-CMIPC-M3C03_G16_s20180431702207_e20180431704580_c20180431705056.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
G = np.sqrt(C.variables['CMI'][:].data) # Band 3 is "green" (0.865 um)

C_file = 'Documents/Delaware/goesR/satdata/OR_ABI-L2-CMIPC-M3C01_G16_s20180431702207_e20180431704580_c20180431705058.nc'  # GOES16 East
C = Dataset(C_file, 'r')
# Load the RGB arrays and apply a gamma correction (square root)
B = np.sqrt(C.variables['CMI'][:].data) # Band 1 is blue (0.47 um)

# "True Green" is some linear interpolation between the three channels
G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

add_seconds = C.variables['t'][0]
DATE = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)


plt.figure(figsize=[10, 8])
plt.imshow(RGB)
plt.title(DATE)



# Satellite height
sat_h = C.variables['goes_imager_projection'].perspective_point_height

# Satellite longitude
sat_lon = C.variables['goes_imager_projection'].longitude_of_projection_origin

# Satellite sweep
sat_sweep = C.variables['goes_imager_projection'].sweep_angle_axis

# The projection x and y coordinates equals
# the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
X = C.variables['x'][:] * sat_h
Y = C.variables['y'][:] * sat_h


# The geostationary projection is perhaps the easiest way, and we don't need to use the Proj object.
# Essentially, we are stretching the image across a map with the same projection and dimensions.
plt.figure(figsize=[10, 8])
m = Basemap(projection='geos', lon_0=sat_lon,
            resolution='i', area_thresh=1000,
            llcrnrx=X.min(),llcrnry=Y.min(),
            urcrnrx=X.max(),urcrnry=Y.max())
m.imshow(np.flipud(RGB)) # Remember, "images" are upside down, so flip up/down
m.drawcoastlines()
m.drawcountries()
m.drawstates()
plt.title(DATE)


# map object with pyproj
p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
# Convert map points to latitude and longitude with the magic provided by Pyproj
XX, YY = np.meshgrid(X, Y)
lons, lats = p(XX, YY, inverse=True)

# Make a new map object for the HRRR model domain map projection
mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
           llcrnrlat=lat0-1.25,llcrnrlon=lon0-1.5,
           urcrnrlat=lat0+1.25,urcrnrlon=lon0+1.5,resolution='h')

xH, yH = mH(lons, lats)


# Create a color tuple for pcolormesh
rgb = RGB[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
rgb = np.minimum(rgb, 1) # Force the maximum possible RGB value to be 1 (the lowest should be 0).
colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]), 3) # flatten array, becuase that's what pcolormesh wants.
colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster?? according to stackoverflow.

# adding this additional line here to see if it helps
# for some reason this helps
colorTuple[colorTuple < 0] = 0
colorTuple[colorTuple > 1] = 1

# Now we can plot the GOES data on the HRRR map domain and projection
plt.figure(figsize=[20, 16])

# The values of R are ignored becuase we plot the color in colorTuple, but pcolormesh still needs its shape.
newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
newmap.set_array(None) # without this line, the linewidth is set to zero, but the RGB colorTuple is ignored. I don't know why.

mH.drawstates()
mH.drawcountries()
mH.drawcoastlines()

plt.title('GOES-16 True Color\n%s' % DATE.strftime('%B %d, %Y %H:%M UTC'))



loc = pyart.io.nexrad_common.get_nexrad_location("KDOX")

lon0 = loc[1] ; lat0 = loc[0]
#use boto to connect to the AWS nexrad holdings directory
s3conn = boto.connect_s3()
bucket = s3conn.get_bucket('noaa-nexrad-level2')
#create a datetime object for the current time in UTC and use the
# year, month, and day to drill down into the NEXRAD directory structure.
now = datetime.utcnow()
date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
        "{:2d}".format(now.day) + '/')
#get the bucket list for the selected date
#Note: this returns a list of all of the radar sites with data for 
# the selected date
ls = bucket.list(prefix=date,delimiter='/')
for key in ls:
    #only pull the data and save the arrays for the site we want
    if site in key.name.split('/')[-2]:
        #set up the path to the NEXRAD files
        path = date + site + '/' + site
        #grab the last file in the file list
        fname = bucket.get_all_keys(prefix=path)[-1]
        #get the file 
        s3key = bucket.get_key(fname)
        #save a temporary file to the local host
        localfile = tempfile.NamedTemporaryFile(delete=False)
        #write the contents of the NEXRAD file to the temporary file
        s3key.get_contents_to_filename(localfile.name)
        #use the read_nexrad_archive function from PyART to read in NEXRAD file
        radar = pyart.io.read_nexrad_archive("Documents/Delaware/goesR/nexraddata/KDOX20180212_171032_V06")
        #get the date and time from the radar file for plot enhancement
        time = radar.time['units'].split(' ')[-1].split('T')
        print(site + ': ' + time[0] + ' at ' + time[1] )

        #set up the plotting grid for the data
        display = pyart.graph.RadarMapDisplay(radar)
        x,y = display._get_x_y(0,True,None)



plt.title('GOES-16 True Color\n%s' % DATE.strftime('%B %d, %Y %H:%M UTC'))

fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(9,9),dpi=100)
#set up a basemap with a lambert conformal projection centered 
# on the radar location, extending 1 degree in the meridional direction
# and 1.5 degrees in the longitudinal in each direction away from the 
# center point.
mH = Basemap(projection='lcc',lon_0=lon0,lat_0=lat0,
           llcrnrlat=lat0-1.4,llcrnrlon=lon0-1.75,
           urcrnrlat=lat0+1.5,urcrnrlon=lon0+1.75,resolution='h')
           
newmap = mH.pcolormesh(xH, yH, R, color=colorTuple, linewidth=0)
newmap.set_array(None) # without this, the linewidth is set to zero, but the RGB color is ignored

#get the plotting grid into lat/lon coordinates
x0,y0 = mH(lon0,lat0)
glons,glats = mH((x0+x*1000.), (y0+y*1000.),inverse=True)
#read in the lowest scan angle reflectivity field in the NEXRAD file 
refl = np.squeeze(radar.get_field(sweep=0,field_name='reflectivity'))
# Mask points with no reflectivity
dBZ = refl
dBZ = np.ma.array(dBZ)
dBZ[dBZ == -10] = np.ma.masked
#set up the plotting parameters (NWSReflectivity colormap, contour levels,
# and colorbar tick labels)
cmap = 'pyart_NWSRef'
levs = np.linspace(0,80,41,endpoint=True)
ticks = np.linspace(0,80,9,endpoint=True)
label = 'Radar Reflectivity Factor ($\mathsf{dBZ}$)'
#define the plot axis to the be axis defined above
ax = axes
#normalize the colormap based on the levels provided above
norm = mpl.colors.BoundaryNorm(levs,256)
cs = mH.pcolormesh(glons,glats,dBZ,norm=norm,cmap=cmap,ax=ax,latlon=True)
#add geographic boundaries and lat/lon labels
mH.drawparallels(np.arange(20,70,0.5),labels=[1,0,0,0],fontsize=12,
                color='k',ax=ax,linewidth=0.001)
mH.drawmeridians(np.arange(-150,-50,1),labels=[0,0,1,0],fontsize=12,
               color='k',ax=ax,linewidth=0.001)
mH.drawcounties(linewidth=0.5,color='k',ax=ax)
mH.drawstates(linewidth=1.5,color='k',ax=ax)
mH.drawcoastlines(linewidth=1.5,color='k',ax=ax)
#mark the radar location with a black dot
mH.scatter(lon0,lat0,marker='o',s=20,color='k',ax=ax,latlon=True)
mH.scatter(-75.7506,39.6780,marker='*',s=20,color='k',ax=ax,latlon=True)
#add the colorbar axes and create the colorbar based on the settings above
cax = fig.add_axes([0.075,0.075,0.85,0.025])
cbar = plt.colorbar(cs,ticks=ticks,norm=norm,cax=cax,orientation='horizontal')
cbar.set_label(label,fontsize=12)
cbar.ax.tick_params(labelsize=11)
#add a title to the figure
fig.text(0.5,0.92, site + ' (0.5$^{\circ}$) Reflectivity\n ' + 
        time[0] + ' at ' + time[1],horizontalalignment='center',fontsize=16)
#display the figure
plt.show()


# close and delete the temporary file holding the radar data
localfile.close()
os.remove(localfile.name)


