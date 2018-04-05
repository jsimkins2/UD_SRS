source /home/james/.bash_profile
source activate goes16

dir=/home/sat_ops/goes_r/cloud_prod/noaa_format/data/
(cd $dir ; ls *C01*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt
(cd $dir ; ls *C02*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt
(cd $dir ; ls *C03*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt
(cd $dir ; ls *C13*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C13_logfile.txt

python /home/sat_ops/goes_r/cloud_prod/noaa_format/conus_goes16_TC.py

# Make the gif from the PNGs we just created above
ls /home/sat_ops/goes_r/cloud_prod/noaa_format/image_conus/ | tail -72 > /home/sat_ops/goes_r/cloud_prod/noaa_format/img_logfile.txt

python /home/sat_ops/goes_r/cloud_prod/noaa_format/gifmaker_goes16.py

scp /home/sat_ops/goes_r/cloud_prod/noaa_format/conus_goes16.gif /var/www/html/imagery/


# Now make a midatlantic version of the code above

python /home/sat_ops/goes_r/cloud_prod/noaa_format/MidAtlantic_goes16TC.py

ls /home/sat_ops/goes_r/cloud_prod/noaa_format/image_midatlantic/ | tail -72 > /home/sat_ops/goes_r/cloud_prod/noaa_format/MAimg_logfile.txt

python /home/sat_ops/goes_r/cloud_prod/noaa_format/gifmaker_MidAtlantic.py

scp /home/sat_ops/goes_r/cloud_prod/noaa_format/MidAtlantic_goes16_TC.gif /var/www/html/imagery/

source deactivate
