# note that the file responsible for clearing out all of the files in image_nxrd_goes is cutoff_epoch.sh in sat_ops/goes16/cloud_proc/noaa_format
source /home/james/.bash_profile
source activate goes16

dir=/home/sat_ops/goes_r/cloud_prod/noaa_format/data/
(cd $dir ; ls *C01*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt
(cd $dir ; ls *C02*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt
(cd $dir ; ls *C03*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt

dir=/home/sat_ops/goes_r/nexrad/data/
(cd $dir ; ls *) | tail -30 > /home/sat_ops/goes_r/nexrad/data_nexrad.txt

python /home/sat_ops/goes_r/nexrad/download_nexrad.py
python /home/sat_ops/goes_r/nexrad/nexrad_goes16.py

ls /home/sat_ops/goes_r/nexrad/image_nxrd_goes/ | tail -10 > /home/sat_ops/goes_r/nexrad/imggoesnxrd_logfile.txt

python /home/sat_ops/goes_r/nexrad/gifmaker_goes16_nexrad.py

source deactivate

scp /home/sat_ops/goes_r/nexrad/goes16_nexrad_kdox.gif /var/www/html/imagery/
