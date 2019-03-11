source /home/james/.bash_profile
source activate goes16

goesdir=/home/sat_ops/goes_r/cloud_prod/noaa_format/data/
(cd $goesdir ; ls *C01*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C01_logfile.txt
(cd $goesdir ; ls *C02*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C02_logfile.txt
(cd $goesdir ; ls *C03*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C03_logfile.txt
(cd $goesdir ; ls *C13*) | tail -30 > /home/sat_ops/goes_r/cloud_prod/noaa_format/C13_logfile.txt

A=$(python /home/sat_ops/goes_r/nexrad/check_site.py)
B="${A: -1}"
if [ "B" == "You have mail in /var/spool/mail/james" ]; then
     A=$(python /home/sat_ops/goes_r/nexrad/check_site.py)
     B="${A: -1}"
fi

if [ "$B" == "1" ]; then
     nexdir=/home/sat_ops/goes_r/nexrad/data_kdox/
     image_dir=/home/sat_ops/goes_r/nexrad/image_kdox_goes
fi

if [ "$B" == "2" ]; then
     nexdir=/home/sat_ops/goes_r/nexrad/data_kdix/
     image_dir=/home/sat_ops/goes_r/nexrad/image_kdix_goes
fi

if [ "$B" == "3" ]; then
     nexdir=/home/sat_ops/goes_r/nexrad/data_klwx/
     image_dir=/home/sat_ops/goes_r/nexrad/image_klwx_goes
fi

ls $nexdir | tail -30 > /home/sat_ops/goes_r/nexrad/data_nexrad.txt

python /home/sat_ops/goes_r/nexrad/download_nexrad.py $B
python /home/sat_ops/goes_r/nexrad/kdox_nexrad_goes16.py $B

ls $image_dir | tail -15 > /home/sat_ops/goes_r/nexrad/imggoesnxrd_logfile.txt

python /home/sat_ops/goes_r/nexrad/gifmaker_nexrad_goes16.py $B

source deactivate

scp /home/sat_ops/goes_r/nexrad/kdox_nexrad_goes16.gif /var/www/html/imagery/
