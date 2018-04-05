
find /home/sat_ops/goes_r/cloud_prod/ -type f -name '3*' -mmin -220 -exec cp {} /home/sat_ops/goes_r/cloud_prod/noaa_format/data/ \;

FILES=/home/sat_ops/goes_r/cloud_prod/noaa_format/data/3*
for f in $FILES
do
  temp1=${f#*.}
  temp2=${temp1%_e*}
  temp3=${temp2::(-3)}
  if [[ ! -e "/home/sat_ops/goes_r/cloud_prod/noaa_format/data/$temp3.nc" ]]; then
    mv $f "/home/sat_ops/goes_r/cloud_prod/noaa_format/data/$temp3.nc"
  else 
    rm -f $f
  fi
done

find /home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR* -type f -mmin +500 -exec rm -f {} \;


find /home/sat_ops/goes_r/cloud_prod/noaa_format/image_conus/* -type f -mmin +500 -exec rm -f {} \;
find /home/sat_ops/goes_r/cloud_prod/noaa_format/image_midatlantic/* -type f -mmin +500 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/image_kdox_goes/* -tiype f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/image_kdix_goes/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/image_klwx_goes/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/data_kdox/* -tiype f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/data_kdix/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/data_klwx/* -type f -mmin +220 -exec rm -f {} \;
