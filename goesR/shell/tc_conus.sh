source /home/james/.bash_profile
source activate goes16

dir=/home/sat_ops/goes_r/cloud_prod/noaa_format/data/
(cd $dir ; ls *C01*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt
(cd $dir ; ls *C02*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt
(cd $dir ; ls *C03*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt

python /home/sat_ops/goes_r/cloud_prod/noaa_format/goes16_conus_tc_png.py

ls /home/sat_ops/goes_r/cloud_prod/noaa_format/image_conus/ | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/img_logfile.txt

python /home/sat_ops/goes_r/cloud_prod/noaa_format/gifmaker_goes16.py

source deactivate

scp /home/sat_ops/goes_r/cloud_prod/noaa_format/conus_goes16.gif /var/www/html/imagery/
