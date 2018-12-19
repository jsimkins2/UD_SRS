
source /home/james/.bash_profile
source activate cartopy35


python /home/sat_ops/goesR/data/sst/scripts/save_raw_sst.py
Rscript /home/sat_ops/goesR/data/sst/scripts/reproject_sst.R
python /home/sat_ops/goesR/data/sst/scripts/save_reprojected_sst.py

#python /home/sat_ops/goesR/radar/klwxtest.py













