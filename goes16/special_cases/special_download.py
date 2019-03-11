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

jday = range(254, 261)
hour = range(0,24)

for j in jday:
    for h in hour:
        # hours has 0 in the filename
        if len(str(h)) == 1:
            hourstr = '0' + str(h)
        
        if len(str(h)) == 2:
            hourstr = str(h)
        
                    
        blobs = bucket.list_blobs(prefix='ABI-L2-MCMIPC/'+ str(2018) + '/' + str(j) + '/' + hourstr + '/')
        print blobs
        results = []
        for i in blobs:
            results.append(i)
            for r in results:
                # parse for the filename we want
                temname = str(r)[58:-1]
                print temname
                # if the file already exists, do NOT download and overwrite it
                # adding the 3223350605. so parsing works downstream, should probably change this later
                if os.path.isfile("/home/sat_ops/goesR/special_cases/florence/data/" + temname.split("_e")[0] + ".nc") == False:
                    # call the individual file we want
                    goesfile= bucket.get_blob('ABI-L2-MCMIPC/'+ str(2018) + '/' + str(j) + '/' + hourstr + '/' + temname)
                    # download said file and keep original naming structure
                    goesfile.download_to_filename("/home/sat_ops/goesR/special_cases/florence/data/" + temname.split("_e")[0] + ".nc") 
                    print "Downloading " + temname
