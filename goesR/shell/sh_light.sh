source /home/james/.bash_profile
source activate goes16

FILES=/home/sat_ops/goes_r/lightning/raw_data/OR*
for f in $FILES
do
  temp1=${f#*data/}
  temp2=${temp1%_e*}
  if [[ ! -e "/home/sat_ops/goes_r/lightning/polished_data/$temp2.nc" ]]; then
    mv $f "/home/sat_ops/goes_r/lightning/polished_data/$temp2.nc"
  fi
done

ls /home/sat_ops/goes_r/lightning/polished_data/ | tail -60 > /home/sat_ops/goes_r/lightning/logfile.txt

python /home/sat_ops/goes_r/lightning/lightning_CONUS.py 

scp /home/sat_ops/goes_r/lightning/lightning_conus.gif /var/www/html/imagery/

source deactivate
