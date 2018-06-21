source /home/james/.bash_profile
source activate goes16

FILES=/home/sat_ops/goes_r/night_scan/raw_data/OR*
for f in $FILES
do
  temp1=${f#*data/}
  temp2=${temp1%_e*}
  if [[ ! -e "/home/sat_ops/goes_r/night_scan/polished_data/$temp2.nc" ]]; then
    mv $f "/home/sat_ops/goes_r/night_scan/polished_data/$temp2.nc"
  fi
done

ls /home/sat_ops/goes_r/night_scan/polished_data/ | tail -100 > /home/sat_ops/goes_r/night_scan/logfile.txt

python /home/sat_ops/goes_r/night_scan/conus_tc_nighttimeIR.py

ls /home/sat_ops/goes_r/night_scan/image_conus/ | tail -70 > /home/sat_ops/goes_r/night_scan/img_logfile.txt

python /home/sat_ops/goes_r/night_scan/gifmaker_ns_conus.py

scp /home/sat_ops/goes_r/night_scan/enhanced_conus_tc.gif /var/www/html/imagery/

source deactivate
~                        
