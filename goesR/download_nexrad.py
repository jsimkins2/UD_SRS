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
import pyart.io
import sys
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)

num=1
num = int(sys.argv[1])

if num==1:
    site = 'KDOX'
    image_dir = 'image_kdox_goes/'
    data_dir = 'data_kdox/'

if num==2:
    site = 'KDIX'
    image_dir = 'image_kdix_goes/'
    data_dir = 'data_kdix/'

if num==3:
    site = 'KLWX'
    image_dir = 'image_klwx_goes/'
    data_dir = 'data_klwx/'

print site
'''
# find the most recent file in the data folder 
import glob
import os
list_of_files = glob.glob('/home/sat_ops/goes_r/nexrad/' + data_dir +'*') # * means all if need specific format then *.csv
if len(list_of_files) > 0:
    latest_file = min(list_of_files, key=os.path.getctime)
    latest_file = latest_file[-16:-3]
else:
    latest_file = "20180404_1319"

latest_time = datetime.strptime(latest_file, '%Y%m%d_%H%M')
now = datetime.utcnow()
kdifftime = now - latest_time

print divmod(kdifftime.days * 86400 + kdifftime.seconds, 60)[0]

seq = range(1,6)
# if we haven't pulled data for a site in over an hour, grab a whole lot of data
if divmod(kdifftime.days * 86400 + kdifftime.seconds, 60)[0] > 30:
'''

seq = range(1, 22)

nxlist = []
for s in seq:
    nxlist.append(s*-1)

now_time = datetime.utcnow()
nex_names = []
print "Downloading"
# get the datetimes for each file
for j in xrange(0, len(nxlist)):
    #get the radar location (this is used to set up the basemap and plotting grid)
    loc = pyart.io.nexrad_common.get_nexrad_location(site)
    lon0 = loc[1] ; lat0 = loc[0]
    #use boto to connect to the AWS nexrad holdings directory
    s3conn = boto.connect_s3()
    bucket = s3conn.get_bucket('noaa-nexrad-level2')
    #create a datetime object for the current time in UTC and use the
    # year, month, and day to drill down into the NEXRAD directory structure.
    date = ("{:4d}".format(now_time.year) + '/' + "{:02d}".format(now_time.month) + '/' +
            "{:02d}".format(now_time.day) + '/')
    ls = bucket.list(prefix=date,delimiter='/')
    for key in ls:
        #only pull the data and save the arrays for the site we want
        if site in key.name.split('/')[-2]:
            #set up the path to the NEXRAD files
            path = date + site + '/' + site
            #grab the last file in the file list
            nex_names = bucket.get_all_keys(prefix=path)[nxlist[j]]
            nex = str(nex_names)[-20:-7]
            print nex
            tem = datetime.strptime(nex, '%Y%m%d_%H%M')
            fname = bucket.get_all_keys(prefix=path)[nxlist[j]]
            #get the file 
            s3key = bucket.get_key(fname)
            #save a temporary file to the local host
            localfile = tempfile.NamedTemporaryFile(delete=False)
            #write the contents of the NEXRAD file to the temporary file
            s3key.get_contents_to_filename(localfile.name)
            #use the read_nexrad_archive function from PyART to read in NEXRAD file
            radar = pyart.io.read_nexrad_archive(localfile.name)
            outfile = '/home/sat_ops/goes_r/nexrad/' + data_dir + nex + '.nc'
            if os.path.isfile(outfile) == False:
                pyart.io.write_cfradial(outfile, radar)
