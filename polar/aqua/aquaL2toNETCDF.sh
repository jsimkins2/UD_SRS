# process the l2 files
ls *L2* > filenames.txt
FILES=`cat filenames.txt`

for f in $FILES
do
    gpt.sh -e aqua.reproject.xml -f netCDF4-CF -t "reprojected_$f" $f
done

ls *reprojected* > filenames.txt 
gpt.sh -e aqua.mosaic.xml -f netCDF4-CF -t mosaiced_fullday.nc `cat filenames.txt`

# now process the RGB files
ls *rgb.hdf* > filenames.txt
gpt.sh -e rgb.aqua.mosaic.xml -f netCDF4-CF -t mosaiced_fullday_rgb.nc `cat filenames.txt`


# now use the collocate tool to combine the files
gpt.sh collocate.xml -Smaster=mosaiced_fullday.nc -Sslave=mosaiced_fullday_rgb.nc -t combined_fullday.nc

