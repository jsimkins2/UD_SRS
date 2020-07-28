source /home/james/.bash_profile
source activate maps37

python3 /home/sat_ops/deos/scripts/deosAG_pointsToGrid.py
python3 /home/sat_ops/deos/scripts/static_AgWx_maps.py

python3 /home/sat_ops/deos/scripts/county_AgWx.py

