cd /data/dineof_in
for file in /home/james/roffs/*nc
do
  filename=`basename $file`
  echo $filename
  masknum=`echo $file | cut -d'_' -f2`
  echo $masknum
  cp script/area.template roffs_dineof.init
  sed -i s/\${filename}/$filename/g roffs_dineof.init
  sed -i s/\${masknum}/$masknum/g roffs_dineof.init

  echo "STARTING AT: " `date`
  /data/dineof_in/dineof-3.0/bin/Linux/dineof-3.0-x64-linux roffs_dineof.init
  echo "FINISHED AT: " `date`
  dest_loc=`echo $filename | sed 's/.*_\(area[0-9]\).*/\1_done/'`
  mv /home/james/roffs/$filename input/${dest_loc}
  #/usr/bin/Rscript /home/sat_ops/dineof/temp/lastFrame.R file=$filename

done


