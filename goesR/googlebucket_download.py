from google.cloud import storage # for some reason this takes forever to do on basin, had trouble installing module
from datetime import datetime, timedelta
import os.path
import os

# open a client 
client = storage.Client()

# ping the bucket
bucket = client.get_bucket('gcp-public-data-goes-16')
# example of pulling from the blob
#blob = bucket.get_blob('ABI-L2-CMIPC/2018/092/19/OR_ABI-L2-CMIPC-M3C02_G16_s20180921912200_e20180921914573_c20180921915082.nc')

# get the current time and julian date so we can figure out what folder we want to be in
now = datetime.utcnow()
tt = now.timetuple()

# julian day has 0s in the filename
if len(str(tt.tm_yday)) == 1:
    jday = '00' + str(tt.tm_yday)

if len(str(tt.tm_yday)) == 2:
    jday = '0' + str(tt.tm_yday) 

if len(str(tt.tm_yday)) == 3:
    jday = str(tt.tm_yday)

# hours has 0 in the filename
if len(str(tt.tm_hour)) == 1:
    hourstr = '0' + str(tt.tm_hour)

if len(str(tt.tm_hour)) == 2:
    hourstr = str(tt.tm_hour)


# get a list of the files in folder
blobs = bucket.list_blobs(prefix='ABI-L2-CMIPC/'+ str(now.year) + '/' + jday + '/' + hourstr + '/')
results = []
for i in blobs:
    results.append(i)
    for j in results:
        # parse for the filename we want
        temname = str(j)[57:-1]
        # if the file already exists, do NOT download and overwrite it
        # adding the 3223350605. so parsing works downstream, should probably change this later
        if os.path.isfile("/home/sat_ops/goes_r/cloud_prod/3223350605." + temname) == False:
            # call the individual file we want
            goesfile= bucket.get_blob('ABI-L2-CMIPC/'+ str(now.year) + '/' + jday + '/' + hourstr + '/' + temname)
            # download said file and keep original naming structure
            goesfile.download_to_filename("/home/sat_ops/goes_r/cloud_prod/3223350605." + temname) 
            print "Downloading " + temname

blobs = bucket.list_blobs(prefix='GLM-L2-LCFA/'+ str(now.year) + '/' + jday + '/' + hourstr + '/')
results = []
for i in blobs:
    results.append(i)
    for j in results:
        # parse for the filename we want
        temname = str(j)[56:-1]
        print temname
        # if the file already exists, do NOT download and overwrite it
        # adding the 3223350605. so parsing works downstream, should probably change this later
        if os.path.isfile("/home/sat_ops/goes_r/lightning/data/" + temname) == False:
            # call the individual file we want
            goesfile= bucket.get_blob('GLM-L2-LCFA/'+ str(now.year) + '/' + jday + '/' + hourstr + '/' + temname)
            # download said file and keep original naming structure
            goesfile.download_to_filename("/home/sat_ops/goes_r/lightning/data/" + temname) 
            print "Downloading " + temname