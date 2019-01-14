# OK NOW WE LOAD IN THE HRRR Dataset
nowdate = datetime.utcnow()
cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/latest.xml')
dataset_name = sorted(cat.datasets.keys())[-1]
dataset = cat.datasets[dataset_name]
ds = dataset.remote_access(service='OPENDAP')
ds = NetCDF4DataStore(ds)
ds = xr.open_dataset(ds)

# parse the temperature at various heights
tempiso = ds.metpy.parse_cf('Temperature_isobaric')
hlats = tempiso['y'][:]
hlons = tempiso['x'][:]
hproj = tempiso.metpy.cartopy_crs
hproj = Proj(hproj.proj4_init)
wgs84=Proj("+init=EPSG:4326")


# Grab actual values
t850 = tempiso[1][2].values
t925 = tempiso[1][3].values
tempsurf = ds.metpy.parse_cf('Temperature_height_above_ground')
tsurf = tempsurf[1][0].values

# grab the 1000 to 500 millibar thickness lines
ht1000 = ds.metpy.parse_cf("Geopotential_height_isobaric")[0][3]
ht500 = ds.metpy.parse_cf("Geopotential_height_isobaric")[0][0]
thick = ht500 - ht1000

# create empty lat and lon arrays so that they are same size so we can transform projection using Pyproj

hrrrlons, hrrrlats = np.meshgrid(hlons, hlats)
hrrrlons, hrrrlats = transform(hproj, wgs84,hrrrlons, hrrrlats)

# trim the data 
hlonrav = hrrrlons.ravel()
hlatrav = hrrrlats.ravel()
lats = hlatrav[(np.abs(hlatrav - boundinglat[0])).argmin():(np.abs(hlatrav - boundinglat[1])).argmin()]
lons = hlatrav[(np.abs(hlonrav - boundinglon[1])).argmin():(np.abs(hlonrav - boundinglon[0])).argmin()]
# 

# trim the data to save space

hrrr_t850 = trim_data(hrrrlats, lons, ma.getdata(t850), boundinglat, boundinglon)
hrrr_t925 = trim_data(hrrrlats, hrrrlons, ma.getdata(t925), boundinglat, boundinglon)
hrrr_tsurf = trim_data(hrrrlats, hrrrlons, ma.getdata(tsurf), boundinglat, boundinglon)
thick = trim_data(lats, lons, ma.getdata(thick), boundinglat, boundinglon)

rav_lats = lats.ravel()
rav_lons = lons.ravel()
rav_t850 = hrrr_t850.ravel()
rav_t925 = hrrr_t925.ravel()
rav_tsurf = hrrr_tsurf.ravel()
rav_thick = thick.ravel()

#Grid Data using scipy interpolate
grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon)
grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat)
glon,glat = np.meshgrid(grid_lons,grid_lats)
grid850= griddata((rav_lons,rav_lats),rav_t850,(glon,glat),method='linear')
grid925 = griddata((rav_lons,rav_lats),rav_t925,(glon,glat),method='linear')
gridsurf = griddata((rav_lons,rav_lats),rav_tsurf,(glon,glat),method='linear')
gridthick = griddata((rav_lons,rav_lats),rav_thick,(glon,glat),method='linear')



# plot up the raw temperatures
fig = plt.figure(figsize=[16,16])
ax = plt.subplot(1,1,1, projection=ccrs.PlateCarree())
#ax.set_extent((min_lon, max_lon, min_lat, max_lat))
im1 = ax.pcolormesh(hlons.values,hlats.values, tsurf,cmap=mpl.cm.jet, vmin=260, vmax=285,transform=tempiso.metpy.cartopy_crs)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.5)
cb1 = plt.colorbar(im1)
# plot up the hrrr_lons trimmed data/
fig = plt.figure(figsize=[16,16])
ax = plt.subplot(1,1,1, projection=ccrs.PlateCarree())
#ax.set_extent((min_lon, max_lon, min_lat, max_lat))
im1 = ax.pcolormesh(hrrrlons,hrrrlats, hrrr_tsurf,cmap=mpl.cm.jet, vmin=260, vmax=285)#,transform=tempiso.metpy.cartopy_crs)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lakes', '50m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                edgecolor='black', facecolor='none',linewidth=1.5))
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.5)
cb1 = plt.colorbar(im1)

# plot the gridded temperatures 