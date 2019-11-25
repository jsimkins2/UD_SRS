import collections
import geopandas as gpd
import cartopy
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import shapely.geometry as sgeom
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable

def add_insetmap(axes_extent, map_extent, state_name, facecolor, edgecolor, geometry):
    # create new axes, set its projection
    use_projection = ccrs.Mercator()     # preserve shape well
    #use_projection = ccrs.PlateCarree()   # large distortion in E-W for Alaska
    geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
    sub_ax = plt.axes(axes_extent, projection=use_projection)  # normal units
    sub_ax.set_extent(map_extent, geodetic)  # map extents

    # add basic land, coastlines of the map
    # you may comment out if you don't need them
    sub_ax.add_feature(cartopy.feature.LAND)
    sub_ax.coastlines()

    sub_ax.set_title(state_name)

    # add map `geometry` here
    sub_ax.add_geometries([usa['geometry'][ind]], ccrs.PlateCarree(),
                          facecolor=facecolor, edgecolor=edgecolor)
    # +++ more features can be added here +++

    # plot box around the map
    extent_box = sgeom.box(map_extent[0], map_extent[2], map_extent[1], map_extent[3])
    sub_ax.add_geometries([extent_box], ccrs.PlateCarree(), color='none', linewidth=0.05)

cmap = mpl.colors.ListedColormap(['white', 'lightyellow', 'pink', 'orange', 'red'])
bounds = [0,1,10,20,30,40]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)    

data = {
    'New Jersey':  30,
    'Rhode Island':   0,
    'Massachusetts':   2,
    'Connecticut':    1,
    'Maryland':   6,
    'New York':    3,
    'Delaware':    26,
    'Florida':     361,
    'Ohio':  5,
    'Pennsylvania':  26,
    'Illinois':    8,
    'California':  4,
    'Virginia':    6,
    'Michigan':    5,
    'Indiana':    7,
    'North Carolina':  25,
    'Georgia':     5,
    'Tennessee':   3,
    'New Hampshire':   0,
    'South Carolina':  5,
    'Louisiana':   1,
    'Kentucky':   2,
    'Wisconsin':  27,
    'Washington':  0,
    'Alabama':     3,
    'Missouri':    3,
    'Texas':   21,
    'West Virginia':   2,
    'Vermont':     0,
    'Minnesota':  49,
    'Mississippi':   3,
    'Iowa':  6,
    'Arkansas':    1,
    'Oklahoma':    16,
    'Arizona':     0,
    'Colorado':    10,
    'Maine':  2,
    'Oregon':  1,
    'Kansas':  10,
    'Utah':  0,
    'Nebraska':    7,
    'Nevada':  1,
    'Idaho':   0,
    'New Mexico':  0,
    'South Dakota':  1,
    'North Dakota':  0,
    'Montana':     1,
    'Wyoming':      1,
    'Hawaii': 0,
    'Alaska': 0,
    'Washington D.C.':1
}

data = dict(sorted(data.items(), key=lambda x: x[0]))

world = gpd.GeoDataFrame.from_file("Downloads/ne_10m_admin_1_states_provinces_lakes/ne_10m_admin_1_states_provinces_lakes.shp")

usa = world[world.iso_a2 == 'US']
usa['coords'] = usa['geometry'].apply(lambda x: x.representative_point().coords[:])
usa['coords'] = [coords[0] for coords in usa['coords']]
usa['data'] = data.values()

can = world[world.iso_a2 == 'CA']
can['coords'] = can['geometry'].apply(lambda x: x.representative_point().coords[:])
can['coords'] = [coords[0] for coords in can['coords']]
candata = dict(zip(list(can['name']),[0,0,1,0,0,0,0,0,0,0,0,0,0]))
can['data'] = candata.values()

# simple plot
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-125, -66.5, 20, 70])
usa.plot(ax=ax, column='data')
can.plot(ax=ax, column='data')


fig = plt.figure(figsize=(16,10), dpi=100)
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.LambertConformal())
ax.set_extent([-125, -74, 20, 55], ccrs.Geodetic())
ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)
ax.set_title('Survey Response', size=20)

for state in usa['name']:
    edgecolor = 'black'

    try:
        # use the name of this state to get pop_density
        state_dens = data[state]
    except:
        state_dens = 0

    # simple scheme to assign color to each state
    if state_dens < 1:
        facecolor = "white"
    elif state_dens >= 1 and state_dens < 10:
        facecolor = "lightyellow"
    elif state_dens >= 10 and state_dens < 20:
        facecolor = "pink"
    elif state_dens >= 20 and state_dens <30:
        facecolor = "orange"
    elif state_dens >= 30:
        facecolor = "red"

    ind = usa[usa['name'] == state].index[0]
    
    if state in ("Alaska", "Hawaii"):
        # print("state.attributes['name']:", state.attributes['name'])

        state_name = state

        # prep map settings
        # experiment with the numbers in both `_extents` for your best results
        if state_name == "Alaska":
            # (1) Alaska
            map_extent = (-178, -135, 46, 73)    # degrees: (lonmin,lonmax,latmin,latmax)
            axes_extent = (0.04, 0.02, 0.29, 0.24) # axes units: 0 to 1, (LLx,LLy,width,height)

        if state_name == "Hawaii":
            # (2) Hawii
            map_extent = (-162, -152, 15, 25)
            axes_extent = (0.27, 0.02, 0.15, 0.12)

        # add inset maps
        add_insetmap(axes_extent, map_extent, state_name, \
                     facecolor, \
                     edgecolor, \
                     [usa['geometry'][ind]])
    # the other (conterminous) states go here
    else:
        # `state.geometry` is the polygon to plot
        ax.add_geometries([usa['geometry'][ind]], ccrs.PlateCarree(),
                          facecolor=facecolor, edgecolor=edgecolor)
for state in can['name']:
    edgecolor = 'black'

    try:
        # use the name of this state to get pop_density
        state_dens = candata[state]
    except:
        state_dens = 0

    # simple scheme to assign color to each state
    if state_dens < 1:
        facecolor = "white"
    elif state_dens >= 1 and state_dens < 10:
        facecolor = "lightyellow"
    elif state_dens >= 10 and state_dens < 20:
        facecolor = "pink"
    elif state_dens >= 20:
        facecolor = "orange"

    ind = can[can['name'] == state].index[0]

    ax.add_geometries([can['geometry'][ind]], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)

sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm._A = []

cb=fig.colorbar(sm,ax=ax, label='Number of Survey Responses',spacing='proportional',pad=-0.01,shrink=0.7)
cb.ax.set_yticklabels(["0","1","10", "20", "30+", ""])
plt.tight_layout()
plt.savefig("Downloads/connors_plot.png", dpi=100)


'''
for idx, row in can.iterrows():
    plt.annotate(s=row['name'], xy=row['coords'],
                 horizontalalignment='center')
