source /home/map_maker/.bash_profile
source activate maps37

# run the map making script
python3 /home/map_maker/deos_maps/deos_geopandas_ares.py

# remove the .tif files that were created during map_making
rm /home/map_maker/deos_maps/temp/*

# copy images over to centaur 
scp /home/map_maker/deos_maps/imagery/* james@128.175.28.202:/var/www/deos/images/home_map/




