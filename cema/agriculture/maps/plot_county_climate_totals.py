import xarray as xr
from datetime import datetime
import numpy as np
import pandas as pd
import calendar
from datetime import datetime, timedelta, date
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
from matplotlib.offsetbox import AnchoredText
import matplotlib.image as image

# declare paths
shapePaths = "/home/james/mapLayers/"
colorPaths = "/home/james/colorramps/"
tiffolder = "/home/sat_ops/deos/static_tifs/"
my_dpi = 100


# cumulative county maps
# open up the county agwx datasets
# list of dataframes
nowdate=datetime.utcnow()

print("at county maps now")
countydf = ['chester_agwx.nc', 'ncc_agwx.nc', 'kentc_agwx.nc', 'sussex_agwx.nc']
doydf = ['chester_agwx_climatology.nc', 'ncc_agwx_climatology.nc', 'kent_agwx_climatology.nc', 'sussex_agwx_climatology.nc']
co_names = ['Chester', 'New Castle', 'Kent', 'Sussex']
nowtime = datetime.utcnow()
ytd = pd.to_datetime(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d")) -  pd.to_datetime(nowtime)
daysback_dict = dict(zip(['YTD'], [int(np.abs(ytd.days))]))#dict(zip(['YTD', '3 Months', '1 Month', '1 Week', '1 Day'], [int(np.abs(ytd.days)), 90, 30, 7, 1]))
doy_len = 366 if calendar.isleap(nowdate.year) == True else 365

# Precipitation Map
print('county precip map')
for db in daysback_dict.keys():
    for co in range(0,len(countydf)):
        print(countydf[co])
        df = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['dailyprecip']
        time_recent = pd.to_datetime(df.time.values[-1])
        df = df.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        print('right before multiplication')
        df.values = df.values * (0.0393701)
        cf = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['dailyprecip']
        cf = cf.isel(dayofyear=slice(0, doy_len))
        cf.values = cf.values * (0.0393701)
        ncep = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['NCEPstageIVPrecip']
        ncep = ncep.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        ncep.values = ncep.values * (0.0393701)
        ncepClim = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['NCEPstageIVPrecip']
        ncepClim.values = ncepClim.values * (0.0393701)
        ncepClim = ncepClim.isel(dayofyear=slice(0, doy_len))
        #Set X range here:
        left = date(nowdate.year, 1, 1)  #Makes it easy to quickly change the range
        right = date(nowdate.year, 12, 31)
        datelist = pd.date_range(str(str(nowdate.year) + "-01-01"), str(str(nowdate.year) + "-12-31")).tolist()
        print('right before plotting')
        fig = plt.figure(figsize=(12,8))
        # Create subplot of 
        plt.subplot(2, 1, 1)
        plt.plot(datelist, np.nancumsum(cf.values), linestyle='dashed', c="lightseagreen",  label="Climatology Cum. Sum")
        plt.plot(df.time.values, np.nancumsum(df.values), c="darkgreen",  label="Observed Cum. Sum")
        plt.plot(df.time.values, df.values, c="darkgreen", label="Observed (in/day)", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else "-"
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + " (in)")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Precipitation')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Precipitation Totals (in)', labelpad = 2)
        plt.legend()
        plt.grid()
        fig.tight_layout()
        #Create subplot of ncep
        plt.subplot(2,1,2)
        plt.plot(datelist, np.nancumsum(ncepClim.values), linestyle='dashed', c="deepskyblue",  label="Climatology")
        plt.plot(ncep.time.values, np.nancumsum(ncep.values), c="blue",  label="Observed Cum. Sum")
        plt.plot(ncep.time.values, ncep.values, c="blue", label="Observed (in/day)", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + " (in)")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - NCEP Stage IV Precipitation')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Precipitaiton Totals (in)', labelpad = 2)
        plt.legend()
        plt.grid()
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 1450, 75 ,zorder=30, alpha=1)
        fig.tight_layout()
        plt.savefig("/var/www/html/imagery/AgWx/county/" + co_names[co] + "_YTD_precipitation.png",bbox_inches='tight',pad_inches = 0.1,dpi=my_dpi*1.3)
        plt.close()

# HDD and CDD 
nowdate=datetime.utcnow()
for db in daysback_dict.keys():
    for co in range(0,len(countydf)):
        df = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['HDD']
        time_recent = pd.to_datetime(df.time.values[-1])
        df = df.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        cf = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['HDD']
        cf = cf.isel(dayofyear=slice(0, doy_len))
        ncep = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['CDD']
        ncep = ncep.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        ncepClim = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['CDD']
        ncepClim = ncepClim.isel(dayofyear=slice(0, doy_len))
        #Set X range here:
        left = date(nowdate.year, 1, 1)  #Makes it easy to quickly change the range
        right = date(nowdate.year, 12, 31)
        datelist = pd.date_range(str(str(nowdate.year) + "-01-01"), str(str(nowdate.year) + "-12-31")).tolist()


        fig = plt.figure(figsize=(12,8))
        # Create subplot of 
        plt.subplot(2, 1, 1)
        plt.plot(datelist, np.nancumsum(cf.values), linestyle='dashed', c="coral",  label="Climatology")
        plt.plot(df.time.values, np.nancumsum(df.values), c="maroon",  label="Observed Cum. Sum")
        plt.plot(df.time.values, df.values, c="maroon", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + "")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Heating Degree Days')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Days', labelpad = 2)
        plt.legend()
        plt.grid()
        fig.tight_layout()
        #Create subplot of ncep
        plt.subplot(2,1,2)
        plt.plot(datelist, np.nancumsum(ncepClim.values), linestyle='dashed', c="deepskyblue",  label="Climatology")
        plt.plot(ncep.time.values, np.nancumsum(ncep.values), c="blue",  label="Observed Cum. Sum")
        plt.plot(ncep.time.values, ncep.values, c="blue", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + "")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Cooling Degree Days')
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Days', labelpad = 1)
        plt.legend()
        plt.grid()
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 1450, 75 ,zorder=30, alpha=1)
        fig.tight_layout()
        plt.savefig("/var/www/html/imagery/AgWx/county/" + co_names[co] + "_YTD_HDD_CDD.png",bbox_inches='tight',pad_inches = 0.1,dpi=my_dpi*1.3)
        plt.close()
        
# GDD and Energy Density
nowdate=datetime.utcnow()
for db in daysback_dict.keys():
    for co in range(0,len(countydf)):
        df = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['GDD']
        time_recent = pd.to_datetime(df.time.values[-1])
        df = df.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        cf = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['GDD']
        cf = cf.isel(dayofyear=slice(0, doy_len))
        ncep = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + countydf[co])['energyDens']
        ncep = ncep.sel(time=slice(datetime.strptime(str(str(nowtime.year) + '-01-01'), "%Y-%m-%d"), time_recent))
        ncepClim = xr.open_dataset('http://basin.ceoe.udel.edu/thredds/dodsC/' + doydf[co])['energyDens']
        ncepClim = ncepClim.isel(dayofyear=slice(0, doy_len))
        #Set X range here:
        left = date(nowdate.year, 1, 1)  #Makes it easy to quickly change the range
        right = date(nowdate.year, 12, 31)
        datelist = pd.date_range(str(str(nowdate.year) + "-01-01"), str(str(nowdate.year) + "-12-31")).tolist()


        fig = plt.figure(figsize=(12,8))
        # Create subplot of 
        plt.subplot(2, 1, 1)
        plt.plot(datelist, np.nancumsum(cf.values), linestyle='dashed', c="lightseagreen",  label="Climatology")
        plt.plot(df.time.values, np.nancumsum(df.values), c="darkgreen",  label="Observed Cum. Sum")
        plt.plot(df.time.values, df.values, c="darkgreen", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round(np.nansum(df.values) - np.nansum(cf.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1)) + "")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Growing Degree Days', pad=5)
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Days', labelpad = 4)
        plt.legend()
        plt.grid()
        fig.tight_layout()
        #Create subplot of ncep
        plt.subplot(2,1,2)
        plt.plot(datelist, np.nancumsum(ncepClim.values), linestyle='dashed', c="deepskyblue",  label="Climatology")
        plt.plot(ncep.time.values, np.nancumsum(ncep.values), c="blue",  label="Observed Cum. Sum")
        plt.plot(ncep.time.values, ncep.values, c="blue", label="Observed", alpha=0.5)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b-01')) #Allows me to format the date into months & days
        plt.gca().xaxis.set_tick_params(rotation = 0)  #puts the x-axis labels on an angle
        plt.gca().set_xbound(left, right)  #changes the range of the x-axis
        plusminus = "+" if np.round(np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values),1) > 0 else ""
        box_text = str("YTD Difference = " + plusminus + str(np.round((np.nansum(ncep.values) - np.nansum(ncepClim.sel(dayofyear=slice(0,nowdate.timetuple().tm_yday)).values))/(1e9),2)) + " (J^9)")
        text_box = AnchoredText(box_text, frameon=True, loc=7, pad=0.5)
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        plt.gca().add_artist(text_box)
        plt.title(str(nowtime.year) + ' ' + co_names[co] + ' County - DEOS Energy Density', pad=5)
        plt.xlabel('Date', labelpad = 5)
        plt.ylabel('Total Energy (J)', labelpad = 2)
        plt.legend()
        plt.grid()
        im1 = image.imread(shapePaths + "deos_logo.png")
        plt.figimage(im1, 1450, 75 ,zorder=30, alpha=1)
        fig.tight_layout()
        plt.savefig("/var/www/html/imagery/AgWx/county/" + co_names[co] + "_YTD_GDD_Energy.png",bbox_inches='tight',pad_inches = 0.1,dpi=my_dpi*1.3)
        plt.close()
