source /home/james/.bash_profile
source activate cartopy35

#python /home/sat_ops/goesR/data/sst/scripts/makeup_sst.py
python /home/sat_ops/goesR/data/noaa_sst/backfill/backfill.py
Rscript /home/sat_ops/goesR/data/noaa_sst/backfill/reproject_sst.R
python /home/sat_ops/goesR/data/noaa_sst/backfill/save_reprojected_sst.py


