source /home/james/.bash_profile
source activate cartopy35

python3 /home/sat_ops/goesR/github/UD_SRS/hrrr/download_hrrr_netcdf.py

#python /home/sat_ops/goesR/github/UD_SRS/nexrad/precipitation_depiction.py
python3 /home/sat_ops/goesR/github/UD_SRS/nexrad/precipitation_depiction_unidata_thredds.py
scp /home/sat_ops/goesR/radar/prectype/*.gif /var/www/html/imagery

find /home/sat_ops/goesR/radar/prectype/precKDOX -type f -mmin +120 -exec rm -f {} \;
find /home/sat_ops/goesR/radar/prectype/precKDIX -type f -mmin +120 -exec rm -f {} \;
find /home/sat_ops/goesR/radar/prectype/precKLWX -type f -mmin +120 -exec rm -f {} \;

