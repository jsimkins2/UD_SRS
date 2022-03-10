# this script will generate true color images 

source /home/james/.bash_profile

ls /home/sat_ops/goesR/radar/refKDOX/ -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/radar/refKDOX/{} /var/www/html/imagery/pivotal/refkdox.png

ls /home/sat_ops/goesR/radar/velKDOX -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/radar/velKDOX/{} /var/www/html/imagery/pivotal/velkdox.png

ls /home/sat_ops/goesR/truecolor/ltng_conus -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/truecolor/ltng_conus/{} /var/www/html/imagery/pivotal/tc_conus.png

ls /home/sat_ops/goesR/truecolor/ltng_mid -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/truecolor/ltng_mid/{} /var/www/html/imagery/pivotal/tc_mid.png

ls /home/sat_ops/goesR/indbands/midatl/imgband10 -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/indbands/midatl/imgband10/{} /var/www/html/imagery/pivotal/midatl_band10.png

ls /home/sat_ops/goesR/indbands/midatl/imgband13 -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/indbands/midatl/imgband13/{} /var/www/html/imagery/pivotal/midatl_band13.png

ls /home/sat_ops/goesR/meso/b13_1 -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/meso/b13_1/{} /var/www/html/imagery/pivotal/meso1.png

ls /home/sat_ops/goesR/meso/b13_2 -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/meso/b13_2/{} /var/www/html/imagery/pivotal/meso2.png

ls /home/sat_ops/goesR/lightning/conus -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/lightning/conus/{} /var/www/html/imagery/pivotal/lightning_conus.png

ls /home/sat_ops/goesR/lightning/midatl -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/lightning/midatl/{} /var/www/html/imagery/pivotal/lightning_midatl.png

ls /home/sat_ops/goesR/fulldisk/atlhurr13 -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/fulldisk/atlhurr13/{} /var/www/html/imagery/pivotal/atlhurr13.png

ls /home/sat_ops/goesR/radar/lrgmid -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/radar/lrgmid/{} /var/www/html/imagery/pivotal/lrgmid.png

ls /home/sat_ops/goesR/radar/tcmid -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/radar/tcmid/{} /var/www/html/imagery/pivotal/tcmid.png

ls /home/sat_ops/goesR/radar/tcconus -tr | tail -n 1 | xargs -I{} scp /home/sat_ops/goesR/radar/tcconus/{} /var/www/html/imagery/pivotal/tcconus.png

wget http://deos.udel.edu/images/home_map/deos_airT.png -P /home/sat_ops/weather_page/
wget http://www.deos.udel.edu/images/current/dttemps.png -P /home/sat_ops/weather_page/
wget http://www.deos.udel.edu/images/current/rain24.png -P /home/sat_ops/weather_page/
wget http://www.deos.udel.edu/images/current/rh.png -P /home/sat_ops/weather_page/
wget http://www.deos.udel.edu/images/current/ws.png -P /home/sat_ops/weather_page/

mv /home/sat_ops/weather_page/*.png /var/www/html/imagery/pivotal/
mv /home/sat_ops/weather_page/*.gif /var/www/html/imagery/pivotal/

mogrify -resize 200 /var/www/html/imagery/pivotal/*.png

scp /home/sat_ops/goesR/github/UD_SRS/weather_page/index.html /var/www/html/weather_page/



