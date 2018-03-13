
find /home/sat_ops/goes_r/cloud_prod/3* -type f -mmin -220 -exec cp {} /home/sat_ops/goes_r/cloud_prod/noaa_format/data/ \;

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

find /home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR* -type f -mmin +220 -exec rm -f {} \;

find /home/sat_ops/goes_r/cloud_prod/noaa_format/image_conus/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/image_nxrd_goes/* -type f -mmin +220 -exec rm -f {} \;
