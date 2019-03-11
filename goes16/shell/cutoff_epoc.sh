
find /home/sat_ops/goes_r/cloud_prod/ -type f -name '3223350605*' -mmin -220 -exec cp {} /home/sat_ops/goes_r/cloud_prod/noaa_format/data/ \;

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

# remove data that's been sitting for too long
find /home/sat_ops/goes_r/cloud_prod/noaa_format/data/OR* -type f -mmin +800 -exec rm -f {} \;
find /home/sat_ops/goes_r/cloud_prod/3* -type f -mmin +800 -exec rm -f {} \;
find /home/sat_ops/goes_r/lightning/raw_data/OR* -type f -mmin +50 -exec rm -f {} \;
find /home/sat_ops/goes_r/lightning/polished_data/OR* -type f -mmin +50 -exec rm -f {} \;
find /home/sat_ops/goes_r/night_scan/raw_data/OR* -type f -mmin +800 -exec rm -f {} \;
find /home/sat_ops/goes_r/night_scan/polished_data/OR* -type f -mmin +800 -exec rm -f {} \;

# remove images that have been sitting too long
find /home/sat_ops/goes_r/cloud_prod/noaa_format/image_conus/* -type f -mmin +800 -exec rm -f {} \;
find /home/sat_ops/goes_r/cloud_prod/noaa_format/image_midatlantic/* -type f -mmin +800 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/image_kdox_goes/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/image_kdix_goes/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/image_klwx_goes/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/data_kdox/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/data_kdix/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/nexrad/data_klwx/* -type f -mmin +220 -exec rm -f {} \;

# remove individual band imagery that has been there too long
find /home/sat_ops/goes_r/ind_bands/conus/band01/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band02/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band03/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band04/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band05/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band06/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band07/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band08/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band09/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band10/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band11/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band12/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band13/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band14/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band15/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/conus/band16/* -type f -mmin +220 -exec rm -f {} \;

find /home/sat_ops/goes_r/ind_bands/midAtl/band01/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band02/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band03/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band04/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band05/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band06/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band07/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band08/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band09/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band10/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band11/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band12/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band13/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band14/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band15/* -type f -mmin +220 -exec rm -f {} \;
find /home/sat_ops/goes_r/ind_bands/midAtl/band16/* -type f -mmin +220 -exec rm -f {} \;
