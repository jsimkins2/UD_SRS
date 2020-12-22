source /home/james/.bash_profile

wget https://www.nhc.noaa.gov/xgtwo/two_atl_2d0.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/cpc/latest/droughtotlk.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/rtma_ru/latest/series_001/sfct.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/rtma_ru/latest/series_001/sfctd_b.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/warnings/nwshaz.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/

wget https://maps8.pivotalweather.com/maps/ndfd/latest/ndfd_sfctmax.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/ndfd/latest/ndfd_48hqpf.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/ndfd/latest/ndfd_48hsnow.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/ndfd/latest/ndfd_24hgust.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/wpc/latest/wpc_excessive_rainfall_day1.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/cpc/latest/1montemp.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/cpc/latest/610prcp.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/cpc/latest/1monprcp.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/
wget https://maps8.pivotalweather.com/maps/spc/spcd1cat.us_ma.png --no-check-certificate -P /home/sat_ops/weather_page/


month=`date '+%m'`;
year=`date '+%Y'`;

if [[ $month -eq 01 ]]
then
  year="$((year-1))"
  month=12
else
  month="$((month-1))"
fi


wget https://iri.columbia.edu/wp-content/uploads/$year/$month/figure1.png -P /home/sat_ops/weather_page/
wget https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.sprd2.gif -P /home/sat_ops/weather_page/
wget https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/pna.sprd2.gif -P /home/sat_ops/weather_page/
wget https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.sprd2.gif -P /home/sat_ops/weather_page/
wget https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.sprd2.gif -P /home/sat_ops/weather_page/
wget https://psl.noaa.gov/map/images/sst/sst.anom.seasonal.gif -P /home/sat_ops/weather_page/
#mv sst.anom.seasonal.gif sst.anom.seasonal.png
wget https://www.cpc.ncep.noaa.gov/products/precip/CWlink/MJO/ensplume_full.gif -P /home/sat_ops/weather_page/
#mv ensplume_full.gif ensplume_full.png
wget https://psl.noaa.gov/map/images/olr/olr.anom.90day.gif -P /home/sat_ops/weather_page/
#mv olr.anom.90day.gif olr.anom.90day.png
wget https://psl.noaa.gov/map/images/rnl/500z_90b.rnl.gif -P /home/sat_ops/weather_page/
#mv 500z_90b.rnl.gif 500z_90b.rnl.png
wget https://psl.noaa.gov/forecasts/reforecast2/analogs/images/deterministic_168to336hrs_latest.png -P /home/sat_ops/weather_page/
wget https://psl.noaa.gov/forecasts/reforecast2/wx_maps/images/z500_mean_f240_nhsm.png -P /home/sat_ops/weather_page/
wget https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso_advisory/figure02.gif -P /home/sat_ops/weather_page/
#mv figure02.gif figure02.png
wget https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso_advisory/figure06.gif -P /home/sat_ops/weather_page/
#mv figure06.gif figure06.png

mv /home/sat_ops/weather_page/*.png /var/www/html/imagery/pivotal/
mv /home/sat_ops/weather_page/*.gif /var/www/html/imagery/pivotal/
  
mogrify -resize 200 /var/www/html/imagery/pivotal/*.png
