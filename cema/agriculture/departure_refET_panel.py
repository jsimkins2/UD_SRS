#!/usr/bin/env python
# coding: utf-8

# In[185]:


#serve deosAgPage.py --port 5000 --allow-websocket-origin=dev-web.demac.udel.edu:5000
#steps to recreate envirionment...
### 1. Create new conda environment
### 2. conda install -c conda-forge bokeh geoviews holoviews
### 3. conda install -c conda-forge hvplot

from cartopy import crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib as mpl
from datetime import date, datetime, timedelta
import pandas as pd
import panel as pn #version 0.9.5
import xarray as xr #version 0.11.3
import hvplot.xarray #version 0.5.2
import hvplot.pandas #version 0.5.2
import geoviews as gv #version 1.7.0
from io import BytesIO
import calendar
import numpy as np

gv.extension('bokeh')

# Declare bounds of the data
bounds=(-76.2,38.3,-74.85, 40.3)

# read in the refET dataset
dsRefET = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/DEOSAG.nc")
dsRefET = dsRefET.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))

climo = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/deos_doy_climatology.nc")
climo = climo.sel(latitude=slice(bounds[3], bounds[1]), longitude=slice(bounds[0],bounds[2]))


# In[ ]:





# In[186]:


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


# In[297]:


startcolor = '#8B4513'
midcolor = '#FFFFFF'
endcolor = '#008000'
own_cmap1 = mpl.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )

# create lists for select widget to pull from
sum_dict = dict(zip(['Heating Degree Days', 'Cooling Degree Days','DEOS Precip',  'Energy Density', 
                     'Reference Evapotranspiration', 'Growing Degree Days'],
                    ['HDD', 'CDD',  'dailyprecip', 'energyDens','refET', 'GDD']))
mean_dict = dict(zip(['Mean Temperature', 'Max Temperature', 'Min Temperature', 'Mean Wind Speed', 'Mean Dew Point',
                      'Mean Relative Humidity', 'Max Relative Humidity', 'Min Relative Humidity', 
                      'Mean Soil Temperature', 'Max Soil Temperature', 'Min Soil Temperature',
                      'Mean Volumetric Water Content','Max Volumetric Water Content', 'Min Volumetric Water Content',
                      'Mean Solar', 'Mean Wind Direction', 'Wind Gust', 'Min Wind Chill'],
                     ['meanTemp', 'maxTemp', 'minTemp', 'meanWS','meanDP','meanRH', 'maxRH', 'minRH',
                      'meanST', 'maxST', 'minST','meanVWC', 'maxVWC', 'minVWC','meanSolar', 'meanWD',
                      'dailyGust', 'dailyMinWC']))

datasets = list(sum_dict.keys()) + list(mean_dict.keys()) + ['NCEP Stage IV Precip', 'NCEP Stage IV Precip - DEOS RefET']
dataset = pn.widgets.Select(name='Dataset', options=datasets, value=datasets[10])
cmap_keys = ['Drought', 'Jet', 'Viridis','Red-Yellow-Green', 'Red-Yellow-Blue', 'Cool-Warm', 'Spectral']
cmap_values = [own_cmap1, 'jet', 'viridis', 'RdYlGn', 'RdYlBu' , 'coolwarm', 'Spectral']
cmap_dict = dict(zip(cmap_keys,cmap_values))


dataset = pn.widgets.Select(name='Dataset', options=datasets, value=datasets[6])
cmap = pn.widgets.Select(name='Color Ramp', options=cmap_keys, value=cmap_keys[5])
start_date = pn.widgets.DatePicker(name='Start Date', value=(date.today() + timedelta(days=-3)))
end_date = pn.widgets.DatePicker(name='End Date', value=(date.today() + timedelta(days=-1)))
#clb_min = pn.widgets.Spinner(name="Colorbar Minimum (mm day-1)", value=0, step=1, start=-5000, end=5000, width=100)
#clb_max = pn.widgets.Spinner(name="Colorbar Maximum (mm day-1)", value=1, step=1, start=-5000, end=5000, width=100)
    
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


### second dataset lists for the departure maps
datasets2 = list(sum_dict.keys()) + list(mean_dict.keys()) + ['NCEP Stage IV Precip']
dataset2 = pn.widgets.Select(name='Dataset', options=datasets2, value=datasets2[6])
cmap2 = pn.widgets.Select(name='Color Ramp', options=cmap_keys, value=cmap_keys[5])


# In[298]:


# only use the @pn.depends if we are making an interactive plot call WITHOUT a button!
#@pn.depends(dataset.param.value, start_date.param.value, end_date.param.value, clb_min.param.value, clb_max.param.value, cmap.param.value)
def make_plot1(dataset, start_date, end_date,cmap):              # clb_min, clb_max, 
     # Load and format the data
    sDate = datetime(start_date.year, start_date.month, start_date.day)
    eDate = datetime(end_date.year, end_date.month, end_date.day)

    quad_title = str(str(dataset) + "    Start Date : " +
                    datetime.strftime(sDate + timedelta(days=-1), "%Y-%m-%d") + " - End Date : " + 
                    datetime.strftime(eDate + timedelta(days=-1), "%Y-%m-%d"))
    if any(dataset in s for s in mean_dict.keys()):
        df = dsRefET[mean_dict[dataset]]
        opLabel = 'Avg ' + df.units 
        df = df.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df = df.mean('time')

    if any(dataset in s for s in sum_dict.keys()):
        df = dsRefET[sum_dict[dataset]]
        opLabel = 'Total ' + df.units
        df = df.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df = df.sum('time')
    if dataset == 'NCEP Stage IV Precip':
        df = dsPrec
        opLabel = 'Total mm'
        df = df.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df = df.sum('time')
        dwnldName = str("ncepStageIV.nc")
    if dataset == 'NCEP Stage IV Precip - DEOS RefET':
        dfref = dsRefET['refET']
        df1 = dfref.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df1 = df1.sum('time')
        df2 = dsPrec.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df2 = df2.sum('time')
        df = df2 - df1.values
        df['Precip - ET'] = df.Precipitation_Flux
        df = df.drop('Precipitation_Flux')
        opLabel = 'Total mm'
    
    x = 'longitude'
    y = 'latitude'
     # create the Altair chart object
    chart = df.hvplot.quadmesh(width=width, height=height, x=x, y=y, cmap=cmap_dict[cmap], 
            project=True, geo=True,title=quad_title,xlim=xlim,ylim=ylim,label=opLabel, #clim=(vmin,vmax)
            rasterize=True, dynamic=False) * shp * shp1
    return chart, df

# create update plot window button
def update(event):
    plotwindow1[1].object = make_plot1(dataset.value, start_date.value, end_date.value,cmap.value)[0] # clb_min, clb_max, 

generate_button = pn.widgets.Button(name='Plot', button_type='primary')
generate_button.on_click(update)

############################################################
def download_cb():
    bout = make_plot1(dataset.value, start_date.value, end_date.value,cmap.value)[1].to_netcdf()
    bio = BytesIO()
    bio.write(bout)
    bio.seek(0)
    return bio

fd = pn.widgets.FileDownload(
    callback=download_cb,
    filename='AgWx.nc')

# set the widget box for the widgets to be placed into
sel_box = pn.WidgetBox(dataset, start_date, end_date, cmap, generate_button,
                       pn.layout.Spacer(height=30),fd)

title1       = '<div style="font-size:50px">CEMA Agriculture Dashboard</div>'
instruction = '<div style="font-size:25px">Select a dataset, set your parameters, and click plot</div>'
oggm_logo   = '<a href="http://cema.udel.edu/"><img src="https://lh3.googleusercontent.com/proxy/WDIKz3hvsUgUMPZJpPgUfaznp5BiT-04YPlehRy2BV2HHYCw9xWRH5RwRD3MVCPmcXp6Ouq-8axaYra-KjwjAidNZ4LC" width=170></a>'
pn_logo     = '<a href="https://panel.pyviz.org"><img src="http://panel.pyviz.org/_static/logo_stacked.png" width=140></a>'
cema_logo = '/home/james/agriculture/shapefolder/cema2logo.png'

header = pn.Row(pn.panel(cema_logo, width=170),  pn.layout.Spacer(width=10), 
                pn.Column(pn.Pane(title1, width=1000), pn.Pane(instruction, width=1000)))
plotwindow1 = pn.Row(sel_box, make_plot1(dataset.value, start_date.value, end_date.value, cmap.value)[0]) #clbmin clbmax


# In[299]:


# only use the @pn.depends if we are making an interactive plot call WITHOUT a button!
#@pn.depends(dataset.param.value, start_date.param.value, end_date.param.value, clb_min.param.value, clb_max.param.value, cmap.param.value)
def make_plot2(dataset2, start_date, end_date,cmap2):              # clb_min, clb_max, 
     # Load and format the data
    sDate = datetime(start_date.year, start_date.month, start_date.day)
    eDate = datetime(end_date.year, end_date.month, end_date.day)
    quad_title = str(str(dataset2) + "    Start Date : " +
                    datetime.strftime(sDate + timedelta(days=-1), "%Y-%m-%d") + " - End Date : " + 
                    datetime.strftime(eDate + timedelta(days=-1), "%Y-%m-%d"))
    # special clauses for climatology file which has a dayofyear index instead of time
    if calendar.isleap(sDate.timetuple().tm_year) == False:
        redim_climo = climo.drop([60], dim='dayofyear')
    else:
        redim_climo = climo
    datetime_list = pd.date_range(start='1/1/' + str(sDate.timetuple().tm_year) , end='12/31/' + str(sDate.timetuple().tm_year))
    redim_climo = redim_climo.assign_coords(time=("time", datetime_list))
    redim_climo = redim_climo.reset_index(['dayofyear'], drop = True)
    for var in redim_climo.data_vars:
        redim_climo[var]= redim_climo[var].rename({'dayofyear': 'time'})

    if sDate.timetuple().tm_year < eDate.timetuple().tm_year:
        y_mult = eDate.timetuple().tm_year - sDate.timetuple().tm_year
        redim_list = []
        for di in range(sDate.timetuple().tm_year, eDate.timetuple().tm_year + 1):
            if calendar.isleap(di) == False:
                redim_list.append(climo.drop([60], dim='dayofyear'))
            else:
                redim_list.append(climo)
        redim_climo = xr.concat(redim_list, dim='dayofyear')
        datetime_list = pd.date_range(start='1/1/' + str(sDate.timetuple().tm_year) , end='12/31/' + str(eDate.timetuple().tm_year))
        redim_climo = redim_climo.assign_coords(time=("time", datetime_list))
        redim_climo = redim_climo.reset_index(['dayofyear'], drop = True)
        for var in redim_climo.data_vars:
            redim_climo[var]= redim_climo[var].rename({'dayofyear': 'time'})
        
    if any(dataset2 in s for s in mean_dict.keys()):
        df2 = dsRefET[mean_dict[dataset2]]
        opLabel = "Avg Diff. (" + df2.units + ")"
        df2 = df2.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df2 = df2.mean('time')
        cf = redim_climo[mean_dict[dataset2]]
        cf = cf.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        cf = cf.mean('time')
        df2 = df2 - cf
        tem_vmin = np.min(df2.values)
        tem_vmax = np.max(df2.values)
        if abs(tem_vmin) >= abs(tem_vmax):
            vmin = tem_vmin
            vmax = (-1)*tem_vmin
        if abs(tem_vmax) > abs(tem_vmin):
            vmax = tem_vmax
            vmin = (-1)*tem_vmax

    if any(dataset2 in s for s in sum_dict.keys()):
        df = dsRefET[sum_dict[dataset2]]
        opLabel = "Total Diff. (" + df.units + ")"
        df = df.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df = df.sum('time')
        cf = redim_climo[sum_dict[dataset2]]
        cf = cf.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        cf = cf.sum('time')
        df2 = df - cf
        tem_vmin = np.min(df2.values)
        tem_vmax = np.max(df2.values)
        if abs(tem_vmin) >= abs(tem_vmax):
            vmin = tem_vmin
            vmax = (-1)*tem_vmin
        if abs(tem_vmax) > abs(tem_vmin):
            vmax = tem_vmax
            vmin = (-1)*tem_vmax
    if dataset2 == 'NCEP Stage IV Precip':
        df = dsPrec
        opLabel = "Total Diff. (mm)"
        df = df.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        df = df.sum('time')
        cf = redim_climo['NCEPstageIVPrecip']
        cf = cf.sel(time=slice(sDate + timedelta(days=-1),eDate + timedelta(days=-1)))
        cf = cf.sum('time')
        df2 = df - cf.values
        dwnldName = str("ncepStageIV.nc")
        tem_vmin = np.min(df2.Precipitation_Flux.values)
        tem_vmax = np.max(df2.Precipitation_Flux.values)
        if abs(tem_vmin) >= abs(tem_vmax):
            vmin = tem_vmin
            vmax = (-1)*tem_vmin
        if abs(tem_vmax) > abs(tem_vmin):
            vmax = tem_vmax
            vmin = (-1)*tem_vmax
        

    x = 'longitude'
    y = 'latitude'

     # create the Altair chart object
    chart2 = df2.hvplot.quadmesh(width=width, height=height, x=x, y=y, cmap=cmap_dict[cmap2], 
            project=True, geo=True,title=quad_title,xlim=xlim,ylim=ylim,label=opLabel, clim=(vmin,vmax),
            rasterize=True, dynamic=False) * shp * shp1
    return chart2, df2



# create update plot window button
def update2(event):
    plotwindow2[1].object = make_plot2(dataset2.value, start_date.value, end_date.value,cmap2.value)[0] # clb_min, clb_max, 

generate_button2 = pn.widgets.Button(name='Plot', button_type='primary')
generate_button2.on_click(update2)

############################################################
def download_cb2():
    bout = make_plot2(dataset2.value, start_date.value, end_date.value,cmap2.value)[1].to_netcdf()
    bio = BytesIO()
    bio.write(bout)
    bio.seek(0)
    return bio

fd2 = pn.widgets.FileDownload(
    callback=download_cb2,
    filename='Departure_AgWx.nc')

# set the widget box for the widgets to be placed into
sel_box2 = pn.WidgetBox(dataset2, start_date, end_date, cmap2, generate_button2,
                       pn.layout.Spacer(height=30),fd2)

title2       = '<div style="font-size:50px">CEMA Agricultural Departure Dashboard</div>'
instruction2 = '<div style="font-size:25px">Select a dataset, set your parameters, and click plot</div>'
oggm_logo   = '<a href="http://cema.udel.edu/"><img src="https://lh3.googleusercontent.com/proxy/WDIKz3hvsUgUMPZJpPgUfaznp5BiT-04YPlehRy2BV2HHYCw9xWRH5RwRD3MVCPmcXp6Ouq-8axaYra-KjwjAidNZ4LC" width=170></a>'
pn_logo     = '<a href="https://panel.pyviz.org"><img src="http://panel.pyviz.org/_static/logo_stacked.png" width=140></a>'
cema_logo = '/home/james/agriculture/shapefolder/cema2logo.png'
header2 = pn.Row(pn.panel(cema_logo, width=170),  pn.layout.Spacer(width=10), 
                pn.Column(pn.Pane(title2, width=1000), pn.Pane(instruction2, width=1000)))

plotwindow2 = pn.Row(sel_box2, make_plot2(dataset2.value, start_date.value, end_date.value, cmap2.value)[0]) #clbmin clbmax


# In[300]:


# now add in another option for the line plots
years_list = list(range(2010, date.today().timetuple().tm_year + 1))
start_year = pn.widgets.Select(name='Start Date', options=years_list, value=years_list[-1])
end_year = pn.widgets.Select(name='Start Date', options=years_list, value=years_list[-1])

county_list = ['Chester County', 'New Castle County', 'Kent County', 'Sussex County']
county = pn.widgets.Select(name='County', options=county_list, value=county_list[1])
datasets3 = list(sum_dict.keys()) + list(mean_dict.keys()) + ['NCEP Stage IV Precip']
dataset3 = pn.widgets.Select(name='Variable', options=datasets3, value=datasets3[6])

def make_plot3(dataset3, start_year, end_year, county):              # clb_min, clb_max, 
     # Load and format the data
    if start_year > end_year:
        start_year = end_year
    sDate = datetime(start_year,1, 1)
    eDate = datetime(end_year, 12, 31)
    quad_title = str(county + " " + str(dataset3) + " (" + str(start_year) + " to " +  str(end_year) + ") vs. Climatology")

    if county == 'Chester County':
        county_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/chester_agwx.nc")
        clim_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/chester_agwx_climatology.nc")
    if county == 'New Castle County':
        county_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/ncc_agwx.nc")
        clim_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/ncc_agwx_climatology.nc")
    if county == "Kent County":
        county_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/kentc_agwx.nc")
        clim_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/kent_agwx_climatology.nc")
    if county == "Sussex County":
        county_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/sussex_agwx.nc")
        clim_df = xr.open_dataset("http://basin.ceoe.udel.edu/thredds/dodsC/sussex_agwx_climatology.nc")
    
    # special clauses for climatology file which has a dayofyear index instead of time
    if sDate.timetuple().tm_year == eDate.timetuple().tm_year:
        if calendar.isleap(sDate.timetuple().tm_year) == False:
            clim_df = clim_df.drop([60], dim='dayofyear')
        else:
            clim_df = clim_df
        datetime_list = pd.date_range(start='1/1/' + str(sDate.timetuple().tm_year) , end='12/31/' + str(sDate.timetuple().tm_year))
        clim_df = clim_df.assign_coords(time=("time", datetime_list))
        clim_df = clim_df.reset_index(['dayofyear'], drop = True)
        for var in clim_df.data_vars:
            clim_df[var]= clim_df[var].rename({'dayofyear': 'time'})
    if sDate.timetuple().tm_year < eDate.timetuple().tm_year:
        y_mult = eDate.timetuple().tm_year - sDate.timetuple().tm_year
        redim_list = []
        for di in range(sDate.timetuple().tm_year, eDate.timetuple().tm_year + 1):
            if calendar.isleap(di) == False:
                redim_list.append(clim_df.drop([60], dim='dayofyear'))
            else:
                redim_list.append(clim_df)
        clim_df = xr.concat(redim_list, dim='dayofyear')
        datetime_list = pd.date_range(start='1/1/' + str(sDate.timetuple().tm_year) , end='12/31/' + str(eDate.timetuple().tm_year))
        clim_df = clim_df.assign_coords(time=("time", datetime_list))
        clim_df = clim_df.reset_index(['dayofyear'], drop = True)
        for var in clim_df.data_vars:
            clim_df[var]= clim_df[var].rename({'dayofyear': 'time'})

                     
    if any(dataset3 in s for s in mean_dict.keys()):      
        df3 = county_df[mean_dict[dataset3]]
        yval = mean_dict[dataset3]
        regLabel = "Avg (" + df3.units + ")"
        df3 = df3.sel(time=slice(sDate,eDate))
        cf3 = clim_df[mean_dict[dataset3]]
        cf3 = cf3.sel(time=slice(sDate,eDate))
        climLabel = "Climatology"
        # create the Altair chart object
        chart3 = cf3.hvplot(height=800, width=1200, x="time", y=yval, label = climLabel, grid=True, title=quad_title, legend = 'left', line_width=3, color='midnightblue')  * df3.hvplot(x="time",
                        y=yval, label = regLabel, legend = 'left', line_width=3, color='seagreen')
        
    if any(dataset3 in s for s in sum_dict.keys()):
        sDate = datetime(start_year,1, 1)
        eDate = datetime(end_year, 12, 31)        
        df3 = county_df[sum_dict[dataset3]]
        yval = sum_dict[dataset3]
        regLabel = "Observed - Cum. Sum (" + df3.units + ")"
        df3 = df3.sel(time=slice(sDate,eDate))
        ddf3 = df3.sel(time=slice(sDate,eDate))
        ddf3 = ddf3.cumsum(skipna=True)
        cf3 = clim_df[sum_dict[dataset3]]
        cf3 = cf3.sel(time=slice(sDate,eDate))
        ccf3 = cf3.sel(time=slice(sDate,eDate))
        ccf3.values = ccf3.values.cumsum()
        climLabel = "Climatology - Cum. Sum (mm)"
        chart3 = ccf3.hvplot(height=800, width=1200, x="time", y=yval, label = climLabel, title=quad_title, grid=True,legend = 'left', line_width=3, color='midnightblue')  * ddf3.hvplot(x="time",
                        y=yval, label = regLabel, legend = 'left', line_width=3, color='seagreen') * cf3.hvplot(x="time", y=yval, label= "Climatology (mm/day)", line_width=2, color='midnightblue') * df3.hvplot(x="time",
                        y=yval, label = "Obsereved (mm/day)", line_width=2, color='seagreen')
        
    if dataset3 == 'NCEP Stage IV Precip':
        df3 = county_df['NCEPstageIVPrecip']
        yval = 'NCEPstageIVPrecip'
        regLabel = "Observed - Cum. Sum (mm)"
        df3 = df3.sel(time=slice(sDate,eDate))
        ddf3 = df3.sel(time=slice(sDate,eDate))
        ddf3.values = ddf3.values.cumsum()
        cf3 = clim_df['NCEPstageIVPrecip']
        cf3 = cf3.sel(time=slice(sDate,eDate))
        ccf3 = cf3.sel(time=slice(sDate,eDate))
        ccf3.values = ccf3.values.cumsum()
        climLabel = "Climatology - Cum. Sum (mm)"
        chart3 = ccf3.hvplot(height=800, width=1200, x="time", y=yval, label = climLabel, title=quad_title, grid=True,legend = 'left', line_width=3, color='midnightblue')  * ddf3.hvplot(x="time",
                        y=yval, label = regLabel, legend = 'left', line_width=3, color='seagreen') * cf3.hvplot(x="time", y=yval, label= "Climatology (mm/day)", line_width=2, color='midnightblue') * df3.hvplot(x="time",
                        y=yval, label = "Obsereved (mm/day)", line_width=2, color='seagreen')
     
    return chart3, df3


##### where i left off 
                        
# need to update labels to add county, and then adjust the stuff below
# create update plot window button
def update3(event):
    plotwindow3[1].object = make_plot3(dataset3.value, start_year.value, end_year.value, county.value)[0] # clb_min, clb_max, 

generate_button3 = pn.widgets.Button(name='Plot', button_type='primary')
generate_button3.on_click(update3)

############################################################
def download_cb3():
    bout = make_plot3(dataset3.value, start_year.value, end_year.value, county.value)[1].to_netcdf()
    bio = BytesIO()
    bio.write(bout)
    bio.seek(0)
    return bio

fd3 = pn.widgets.FileDownload(
    callback=download_cb3,
    filename='County_AgWx.nc')

# set the widget box for the widgets to be placed into
sel_box3 = pn.WidgetBox(dataset3, start_year, end_year, county, generate_button3,
                       pn.layout.Spacer(height=30),fd3)

title3       = '<div style="font-size:50px">CEMA County Averaged Agriculture Dashboard</div>'
instruction3 = '<div style="font-size:25px">Select a dataset, set your parameters, and click plot</div>'
oggm_logo   = '<a href="http://cema.udel.edu/"><img src="https://lh3.googleusercontent.com/proxy/WDIKz3hvsUgUMPZJpPgUfaznp5BiT-04YPlehRy2BV2HHYCw9xWRH5RwRD3MVCPmcXp6Ouq-8axaYra-KjwjAidNZ4LC" width=170></a>'
pn_logo     = '<a href="https://panel.pyviz.org"><img src="http://panel.pyviz.org/_static/logo_stacked.png" width=140></a>'
cema_logo = '/home/james/agriculture/shapefolder/cema2logo.png'
header3 = pn.Row(pn.panel(cema_logo, width=170),  pn.layout.Spacer(width=10), 
                pn.Column(pn.Pane(title3, width=1000), pn.Pane(instruction3, width=1000)))

plotwindow3 = pn.Row(sel_box3, make_plot3(dataset3.value, start_year.value, end_year.value, county.value)[0]) #clbmin clbmax


# In[301]:


dashboard = pn.Column(header,pn.layout.Spacer(height=50), plotwindow1, pn.layout.Spacer(height=50),header2,
                      pn.layout.Spacer(height=50), plotwindow2,
                      pn.layout.Spacer(height=50),header3, pn.layout.Spacer(height=50), plotwindow3)
dashboard.servable()

