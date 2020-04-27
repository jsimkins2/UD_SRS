# panel serve refET_panel.py --port 8081 --allow-websocket-origin=agmap.cema.udel.edu:8081

from cartopy import crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib as mpl
from datetime import date, datetime, timedelta
import panel as pn #version 0.9.5
import xarray as xr #version 0.11.3
import hvplot.xarray #version 0.5.2
import hvplot.pandas #version 0.5.2
import geoviews as gv #version 1.7.0
gv.extension('bokeh')

# Declare bounds of the data
bounds=(-76.2,38.3,-74.85, 40.3)

# read in the refET dataset
dsRefET = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSRefET.nc")
dsRefET = dsRefET.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))


# In[2]:


# read in the ncep stage IV precip dataset
dsPrec = xr.open_dataset("http://thredds.demac.udel.edu/thredds/dodsC/NCEPIVQC.nc")
dsPrec = dsPrec.sel(lat=slice(bounds[1], bounds[3]), 
                    lon=slice(bounds[0],bounds[2]), 
                    time=slice(datetime.strptime("2014-01-01", "%Y-%m-%d"),
                              date.today()))

dsPrec = dsPrec.reindex(lat=list(reversed(dsPrec.lat)))
dsPrec = dsPrec.rename(name_dict= {'lat' : 'latitude'})
dsPrec = dsPrec.rename(name_dict= {'lon' : 'longitude'})
dsPrec = dsPrec.drop('crs')


# In[3]:


# define own colorbar
startcolor = '#8B4513'
midcolor = '#FFFFFF'
endcolor = '#008000'
own_cmap1 = mpl.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )

# create lists for select widget to pull from
datasets = ['RefET', 'NCEP Stage IV Precip', 'NCEP Stage IV Precip - RefET']
cmap_keys = ['Drought', 'Jet', 'Viridis','Red-Yellow-Green', 'Red-Yellow-Blue', 'Cool-Warm', 'Spectral']
cmap_values = [own_cmap1, 'jet', 'viridis', 'RdYlGn', 'RdYlBu' , 'coolwarm', 'Spectral']
cmap_dict = dict(zip(cmap_keys,cmap_values))

# define widgets
dataset = pn.widgets.Select(name='Dataset', options=datasets, value=datasets[0])
cmap = pn.widgets.Select(name='Color Ramp', options=cmap_keys, value=cmap_keys[0])
start_date = pn.widgets.DatePicker(name='Start Date', value=(date.today() + timedelta(days=-3)))
end_date = pn.widgets.DatePicker(name='End Date', value=(date.today() + timedelta(days=-1)))
#clb_min = pn.widgets.Spinner(name="Colorbar Minimum (mm day-1)", value=-5000, step=1, start=-5000, end=5000, width=100)
#clb_max = pn.widgets.Spinner(name="Colorbar Maximum (mm day-1)", value=-4999, step=1, start=-5000, end=5000, width=100)
    
# information for quadmesh plot
ylim=(38.3, 40.3)
xlim=(-76.2,-74.85)
width=800
height=1000

# load in the shapefiles
deos = shpreader.Reader("/home/james/agriculture/shapefolder/cb_2018_us_state_500k/" + 'cb_2018_us_state_500kclipped.shp')
coastline = gv.Shape.from_records(deos.records())
shp = gv.Polygons(coastline).opts('Polygons', line_color='black',fill_color=None)

deos = shpreader.Reader("/home/james/agriculture/shapefolder/mapLayers/" + 'deoscounties.shp')
coastline = gv.Shape.from_records(deos.records())
shp1 = gv.Polygons(coastline).opts('Polygons', line_color='black',fill_color=None)


# In[4]:


# only use the @pn.depends if we are making an interactive plot call WITHOUT a button!
#@pn.depends(dataset.param.value, start_date.param.value, end_date.param.value, clb_min.param.value, clb_max.param.value, cmap.param.value)
def make_plot(dataset, start_date, end_date,cmap):              # clb_min, clb_max, 
     # Load and format the data
    sDate = datetime(start_date.year, start_date.month, start_date.day)
    eDate = datetime(end_date.year, end_date.month, end_date.day)

    quad_title = str(str(dataset) + "    Start Date : " +
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
    if dataset == 'NCEP Stage IV Precip - RefET':
        df1 = dsRefET.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df1 = df1.sum('time')
        df2 = dsPrec.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df2 = df2.sum('time')
        df = df2 - df1.refET.values
        df['Precip - ET'] = df.Precipitation_Flux
        df = df.drop('Precipitation_Flux')


    x = 'longitude'
    y = 'latitude'
     # create the Altair chart object
    chart = df.hvplot.quadmesh(width=width, height=height, x=x, y=y, cmap=cmap_dict[cmap], 
            project=True, geo=True,title=quad_title,xlim=xlim,ylim=ylim,label="mm", #clim=(vmin,vmax)
            rasterize=True, dynamic=False) * shp * shp1
    return chart


# In[5]:


def update(event):
    plotwindow[1].object = make_plot(dataset.value, start_date.value, end_date.value,cmap.value) # clb_min, clb_max, 

generate_button = pn.widgets.Button(name='Plot', button_type='primary')
generate_button.on_click(update)

sel_box = pn.WidgetBox(dataset, start_date, end_date, cmap, generate_button) # clb_min, clb_max, 


# In[6]:


## add in the headers
title       = '<div style="font-size:50px">CEMA Agriculture Dashboard</div>'
instruction = '<div style="font-size:25px">Select a dataset, set your parameters, and click plot</div>'
oggm_logo   = '<a href="http://cema.udel.edu/"><img src="https://lh3.googleusercontent.com/proxy/WDIKz3hvsUgUMPZJpPgUfaznp5BiT-04YPlehRy2BV2HHYCw9xWRH5RwRD3MVCPmcXp6Ouq-8axaYra-KjwjAidNZ4LC" width=170></a>'
pn_logo     = '<a href="https://panel.pyviz.org"><img src="http://panel.pyviz.org/_static/logo_stacked.png" width=140></a>'
cema_logo = '/home/james/agriculture/shapefolder/cema2logo.png'

header = pn.Row(pn.panel(cema_logo, width=170),  pn.layout.Spacer(width=10), 
                pn.Column(pn.Pane(title, width=1000), pn.Pane(instruction, width=1000)))


# In[7]:


plotwindow = pn.Row(sel_box, make_plot(dataset.value, start_date.value, end_date.value, cmap.value)) #clbmin clbmax
dashboard = pn.Column(header,pn.layout.Spacer(height=50), plotwindow)
dashboard.servable()


# In[ ]:




