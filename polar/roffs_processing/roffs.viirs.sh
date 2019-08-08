#!/bin/bash
# this is really for level 1 viirs files to be processed to geotiff files
/opt/ocssw/scripts/getanc.py -r SVM01.h5 
raw_name=`ls *.*.*.*`
year=`echo $raw_name | sed -n 's/\([0-9]\{4\}\).*/\1/p'`
month=`echo $raw_name | sed -n 's/[0-9]\{4\}.\([0-9]\{2\}\).*/\1/p'`
day=`echo $raw_name | sed -n 's/[0-9]\{4\}.[0-9]\{2\}\([0-9]\{2\}\).*/\1/p'`
ptime=`echo $raw_name | sed -n 's/[0-9]\{4\}.[0-9]\{4\}.\([0-9]\{6\}\).*/\1/p'`
jday=`date -d "$year-$month-$day" +%j`
l2filename=V${year}${jday}${ptime}.L2_1KM

echo "l2gen ifile=SVM01.h5 geofile=GMTCO.h5 ofile=${l2filename} l2prod1="sst qual_sst l2_flags chl_oc3 a_410_qaa a_443_qaa a_486_qaa a_551_qaa a_671_qaa bb_551_qaa aph_443_qaa adg_410_qaa c_551_qaa Rrs_410 Rrs_443 Rrs_486 Rrs_551 Rrs_671 Rrs_745 Rrs_862 ndvi evi pic poc class_34k_w_owmc" resolution=1000 ctl_pt_incr=1 ctl_pt_incr=1 proc_ocean=1 proc_sst=1 proc_land=0 atmocor=1 maskcloud=1 maskland=1 maskhilt=0 maskstlight=0 par=SVM01.h5.anc"

l2gen ifile=SVM01.h5 geofile=GMTCO.h5 ofile=${l2filename} l2prod1="sst qual_sst l2_flags chl_oc3 a_410_qaa a_443_qaa a_486_qaa a_551_qaa a_671_qaa bb_551_qaa aph_443_qaa adg_410_qaa c_551_qaa Rrs_410 Rrs_443 Rrs_486 Rrs_551 Rrs_671 Rrs_745 Rrs_862 ndvi evi pic poc class_34k_w_owmc" resolution=1000 ctl_pt_incr=1 ctl_pt_incr=1 proc_ocean=1 proc_sst=1 proc_land=0 atmocor=1 maskcloud=1 maskland=1 maskhilt=0 maskstlight=0 sline=3 subsamp=8 par=SVM01.h5.anc

if [[ $? != 0 ]]; then
  l2gen ifile=SVM01.h5 geofile=GMTCO.h5 ofile=${l2filename} l2prod1="sst qual_sst l2_flags chl_oc3 a_410_qaa a_443_qaa a_486_qaa a_551_qaa a_671_qaa bb_551_qaa aph_443_qaa adg_410_qaa c_551_qaa Rrs_410 Rrs_443 Rrs_486 Rrs_551 Rrs_671 Rrs_745 Rrs_862 ndvi evi pic poc class_34k_w_owmc" resolution=1000 ctl_pt_incr=1 ctl_pt_incr=1 proc_ocean=1 proc_sst=1 proc_land=0 atmocor=1 maskcloud=1 maskland=1 maskhilt=0 maskstlight=0 par=SVM01.h5.anc calfile=
fi


# ao1 cover area
gpt.sh -e ao1.reproject.xml -f netCDF4-CF -t ${l2filename}.ao1.r.nc ${l2filename}
gpt.sh -e ao1.mosaic.xml -f netCDF4-CF -t ${l2filename}.ao1.mr.nc ${l2filename}.ao1.r.nc
gpt.sh -e sst.write.ao1.xml -Ssource=${l2filename}.ao1.mr.nc
gpt.sh -e chl.write.ao1.xml -Ssource=${l2filename}.ao1.mr.nc
convert sst.ao1.png sst.ao1.${l2filename}.tiff
convert chl.ao1.png chl.ao1.${l2filename}.tiff
convert sst.ao1.${l2filename}.tiff -set colorspace Gray -separate -average sst.ao1.${l2filename}.tiff
convert chl.ao1.${l2filename}.tiff -set colorspace Gray -separate -average chl.ao1.${l2filename}.tiff

# ao2 cover area
gpt.sh -e ao2.reproject.xml -f netCDF4-CF -t ${l2filename}.ao2.r.nc ${l2filename}
gpt.sh -e ao2.mosaic.xml -f netCDF4-CF -t ${l2filename}.ao2.mr.nc ${l2filename}.ao2.r.nc
gpt.sh -e sst.write.ao2.xml -Ssource=${l2filename}.ao2.mr.nc
gpt.sh -e chl.write.ao2.xml -Ssource=${l2filename}.ao2.mr.nc
convert sst.ao2.png sst.ao2.${l2filename}.tiff
convert chl.ao2.png chl.ao2.${l2filename}.tiff
convert sst.ao2.${l2filename}.tiff -set colorspace Gray -separate -average sst.ao2.${l2filename}.tiff
convert chl.ao2.${l2filename}.tiff -set colorspace Gray -separate -average chl.ao2.${l2filename}.tiff

# ao3 cover area
gpt.sh -e ao3.reproject.xml -f netCDF4-CF -t ${l2filename}.ao3.r.nc ${l2filename}
gpt.sh -e ao3.mosaic.xml -f netCDF4-CF -t ${l2filename}.ao3.mr.nc ${l2filename}.ao3.r.nc
gpt.sh -e sst.write.ao3.xml -Ssource=${l2filename}.ao3.mr.nc
gpt.sh -e chl.write.ao3.xml -Ssource=${l2filename}.ao3.mr.nc
convert sst.ao3.png sst.ao3.${l2filename}.tiff
convert chl.ao3.png chl.ao3.${l2filename}.tiff
convert sst.ao3.${l2filename}.tiff -set colorspace Gray -separate -average sst.ao3.${l2filename}.tiff
convert chl.ao3.${l2filename}.tiff -set colorspace Gray -separate -average chl.ao3.${l2filename}.tiff
