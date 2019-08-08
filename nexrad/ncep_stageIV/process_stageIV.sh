source /home/james/.bash_profile
Rscript /home/james/ncep_stageIV/scripts/ncepStageIVDataConversion.R

find /home/james/ncep_stageIV/unprojected/* -type f -mmin +21600 -exec rm -f {} \;
find /home/james/ncep_stageIV/notNetcdf/* -type f -mmin +21600 -exec rm -f {} \;




