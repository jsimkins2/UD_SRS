
import panel as pn
import pandas as pd
from datetime import date, datetime, timedelta
import xarray as xr
import hvplot.xarray
import hvplot.pandas
import holoviews as hv
from cartopy import crs as ccrs
import cartopy.io.shapereader as shpreader
from shapely.ops import unary_union
import geoviews as gv
import geoviews.feature as gf
from geoviews import opts
import matplotlib as mpl
from shapely.geometry import box, MultiLineString
gv.extension('bokeh')


# Declare bounds of the data
bounds=(-76.2,38.3,-74.85, 40.3)

# read in the refET dataset
dsRefET = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSRefET.nc")
dsRefET = dsRefET.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))

# read in the ncep stage IV precip dataset
dsPrec = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/ncep_stage_iv.nc")
dsPrec = dsPrec.sel(lat=slice(bounds[1], bounds[3]), 
                    lon=slice(bounds[0],bounds[2]), 
                    time=slice(datetime.strptime("2014-01-01", "%Y-%m-%d"),
                              date.today()))

dsPrec = dsPrec.reindex(lat=list(reversed(dsPrec.lat)))
dsPrec = dsPrec.rename(name_dict= {'lat' : 'latitude'})
dsPrec = dsPrec.rename(name_dict= {'lon' : 'longitude'})


title='DEOS Reference Evapotranspiration'
subtitle='2014-2020'

datasets = ['RefET', 'NCEP Stage IV Precip', 'RefET - NCEP Stage IV Precip']
# this creates the dropdown widget
dataset_val = pn.widgets.Select(name='Dataset', options=datasets, value=datasets[2])
start_date_val = pn.widgets.DatePicker(name='Start Date', value=(date.today() + timedelta(days=-8)))
end_date_val = pn.widgets.DatePicker(name='End Date', value=(date.today() + timedelta(days=-1)))
clb_min_val = pn.widgets.Spinner(name="Colorbar Minimum (mm day-1)", value=0, step=1, start=-5000, end=5000, width=100)
clb_max_val = pn.widgets.Spinner(name="Colorbar Maximum (mm day-1)", value=1, step=1, start=-5000, end=5000, width=100)


startcolor = '#8B4513'
midcolor = '#FFFFFF'
endcolor = '#008000'
own_cmap1 = mpl.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )

ylim=(38.3, 40.3)
xlim=(-76.2,-74.85)
width=800
height=1000

sDate = date.today() + timedelta(days=-8)
eDate = date.today() + timedelta(days=-1)
#extents = (xlim[0], xlim[1], ylim[0], ylim[1])
#states = gv.feature.states.data.with_scale('10m')
#stateLines = gv.Shape(unary_union(list(states.intersecting_geometries(extents))))
#shapePaths = "/Users/James/Downloads/cb_2018_us_county_500k/"
#shpfile = shapePaths + 'cb_2018_us_county_500k.shp'
#bounds=(-76.5,38.3,-74.85, 40.3)
#bbox = box(*bounds)
#reader = shpreader.Reader(shpfile)
#coastline = gv.Shape(unary_union([geom.geometry.intersection(bbox) for geom in reader.records()]))
#shp = gv.Polygons(coastline).opts('Polygons', line_color='black',fill_color=None)
deos = shpreader.Reader("/var/www/html/james/ag_page/mapLayers/" + 'cb_2018_us_state_500kclipped.shp')
coastline = gv.Shape.from_records(deos.records())
shp = gv.Polygons(coastline).opts('Polygons', line_color='black',fill_color=None)

deos = shpreader.Reader("/var/www/html/james/ag_page/mapLayers/" + 'deoscounties.shp')
coastline = gv.Shape.from_records(deos.records())
shp1 = gv.Polygons(coastline).opts('Polygons', line_color='black',fill_color=None)


# In[9]:


@pn.depends(dataset_val.param.value, start_date_val.param.value, end_date_val.param.value,
            clb_min_val.param.value, clb_max_val.param.value)
def make_plot(dataset, start_date, end_date, clb_min, clb_max):
     # Load and format the data
    sDate = datetime(start_date.year, start_date.month, start_date.day)
    edate = datetime(end_date.year, end_date.month, end_date.day)
    
    quad_title = str("DEOS " + str(dataset) + "    Start Date : " +
                    datetime.strftime(sDate + timedelta(days=-1), "%Y-%m-%d") + " - End Date : " + 
                    datetime.strftime(eDate + timedelta(days=-1), "%Y-%m-%d"))
    if dataset == 'RefET':
        df = dsRefET
        df = df.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df = df.sum('time')
    if dataset == 'NCEP Stage IV Precip':
        df = dsPrec
        df = df.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df = df.sum('time')
    if dataset == 'RefET - NCEP Stage IV Precip':
        df1 = dsRefET.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df1 = df1.sum('time')
        df2 = dsPrec.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df2 = df2.sum('time')
        df = df1 - df2.Precipitation_Flux.values


     # create filter mask for the dataframe
    vmin=clb_min
    vmax=clb_max
    x = 'longitude'
    y = 'latitude'
     # create the Altair chart object
    chart = df.hvplot.quadmesh(width=width, height=height, x=x, y=y, cmap=own_cmap1, 
             clim=(vmin,vmax),project=True, geo=True,title=quad_title,xlim=xlim,ylim=ylim,label="mm",
             rasterize=True, dynamic=False) * shp * shp1
    return chart


# In[10]:


def update(event):
    dashboard[1].object = make_plot(dataset_val.param.value, start_date_val.param.value, end_date_val.param.value,
            clb_min_val.param.value, clb_max_val.param.value)

generate_button = pn.widgets.Button(name='Plot', button_type='primary')
generate_button.on_click(update)

sel_box = pn.WidgetBox(dataset_val, start_date_val, end_date_val, clb_min_val, clb_max_val,
                                          generate_button)


# In[11]:


## add in the headers
title       = '<div style="font-size:50px">CEMA Agriculture Dashboard</div>'
instruction = '<div style="font-size:25px">Select a dataset, set your parameters, and click plot</div>'
oggm_logo   = '<a href="http://cema.udel.edu/"><img src="https://lh3.googleusercontent.com/proxy/770dmeH3v2ImKDM29UjbEj3vAXe2V8d4fI6oaZDpDegE4TXBaYM4W0HCsAFAZdKnEHyTxxkVw0bZIttIOpCmpJqVpyW2" width=170></a>'
pn_logo     = '<a href="https://panel.pyviz.org"><img src="http://panel.pyviz.org/_static/logo_stacked.png" width=140></a>'

header = pn.Row(pn.Pane(oggm_logo),  pn.layout.Spacer(width=10), 
                pn.Column(pn.Pane(title, width=1000), pn.Pane(instruction, width=1000)))


# In[12]:


dashboard = pn.Column(header,pn.layout.Spacer(height=50), pn.Row(sel_box, make_plot))


# In[13]:


dashboard.servable()


# In[ ]:




