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
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)

seq = range(1, 5)
nxlist = []
for s in seq:
    nxlist.append(s*-1)

now_time = datetime.utcnow()
nex_names = []
# get the datetimes for each file
for j in xrange(0, len(nxlist)):
    site = 'KDOX'
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
            outfile = '/home/sat_ops/goes_r/nexrad/data/kdox_' + nex
            if os.path.isfile(outfile + ".nc") == False:
                pyart.io.write_cfradial(outfile, radar)
