# check which nexrad location we can use
# script designed for basin.ceoe.udel.edu
# James Simkins
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
from pyproj import Proj     
from scipy import ndimage
from scipy import stats
import os.path
import pyart
import boto
import os
import tempfile
import matplotlib as mpl
from dateutil import tz
import time
from time import mktime
import matplotlib.image as image

#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)

# first things first, organize and only use radar that's actually updating
site = 'KDOX'
image_dir = 'image_kdox_goes/'
data_dir = 'data_kdox/'
# first we need to check if KDOX is up and running
#get the radar location (this is used to set up the basemap and plotting grid)
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]
#use boto to connect to the AWS nexrad holdings directory
s3conn = boto.connect_s3()
bucket = s3conn.get_bucket('noaa-nexrad-level2')
#create a datetime object for the current time in UTC and use the
# year, month, and day to drill down into the NEXRAD directory structure.
now = datetime.utcnow()
date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
        "{:02d}".format(now.day) + '/')
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
        fname = bucket.get_all_keys(prefix=path)[-5]
        #get the file 
        s3key = bucket.get_key(fname)
        #save a temporary file to the local host
        localfile = tempfile.NamedTemporaryFile(delete=False)
        #write the contents of the NEXRAD file to the temporary file
        s3key.get_contents_to_filename(localfile.name)
        #use the read_nexrad_archive function from PyART to read in NEXRAD file
        radar = pyart.io.read_nexrad_archive(localfile.name)
        #get the date and time from the radar file for plot enhancement
        ktime = radar.time['units'].split(' ')[-1].split('T')
        checktime = str(ktime[0] + ' ' + ktime[1][:-1])
        checktime = datetime.strptime(checktime, '%Y-%m-%d %H:%M:%S')
        kdifftime = now - checktime
        
        if divmod(kdifftime.days * 86400 + kdifftime.seconds, 60)[0] > 33:
            site = 'KDIX'
            image_dir = 'image_kdix_goes/'
            data_dir = 'data_kdix/'
            print "KDOX is down, using KDIX now"

print site
loc = pyart.io.nexrad_common.get_nexrad_location(site)
lon0 = loc[1] ; lat0 = loc[0]
# if KDIX is down, well we gotta use something I guess 
s3conn = boto.connect_s3()
bucket = s3conn.get_bucket('noaa-nexrad-level2')
#create a datetime object for the current time in UTC and use the
# year, month, and day to drill down into the NEXRAD directory structure.
now = datetime.utcnow()
date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
        "{:02d}".format(now.day) + '/')
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
        radar = pyart.io.read_nexrad_archive(localfile.name)
        #get the date and time from the radar file for plot enhancement
        ktime = radar.time['units'].split(' ')[-1].split('T')
        checktime = str(ktime[0] + ' ' + ktime[1][:-1])
        checktime = datetime.strptime(checktime, '%Y-%m-%d %H:%M:%S')
        kdifftime = now - checktime
        
        if divmod(kdifftime.days * 86400 + kdifftime.seconds, 60)[0] > 21:
            site = 'KLWX'
            image_dir = 'image_klwx_goes/'
            data_dir = 'data_kdox/'
            print "KDIX is down, using KLWX now"
            
def sitecheck():
    if site == 'KDOX':
        return(1)
    if site == 'KDIX':
        return(2)
    if site == 'KLWX':
        return(3)

x = sitecheck()
print x