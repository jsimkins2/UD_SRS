# note that the file responsible for clearing out all of the files in image_nxrd_goes is cutoff_epoch.sh in sat_ops/goes16/cloud_proc/noaa_format
source /home/james/.bash_profile
source activate goes16
python /home/sat_ops/goes_r/nexrad/kdox_radar_loop.py

scp /home/sat_ops/goes_r/nexrad/kdox_radar.gif /var/www/html/imagery/

