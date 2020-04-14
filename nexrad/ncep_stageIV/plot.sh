source /home/james/.bash_profile
python /home/james/ncep_stageIV/scripts/ncep_rolling24_plot.py

scp /home/james/ncep_stageIV/*.png james@basin.ceoe.udel.edu:/var/www/html/imagery/
