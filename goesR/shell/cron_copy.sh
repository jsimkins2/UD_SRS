SHELL=/bin/sh
PATH=/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:/opt/anaconda2/envs/goes16/

40 * * * * sh /home/james/anomalies/aqua.anomalies.sh
*/5 * * * * sh /home/sat_ops/goes_r/cloud_prod/noaa_format/cutoff_epoc.sh
*/15 * * * * sh /home/sat_ops/goes_r/cloud_prod/noaa_format/tc_conus.sh
*/6 * * * * sh /home/sat_ops/goes_r/nexrad/kdox.sh
*/12 * * * * sh /home/sat_ops/goes_r/nexrad/nxrd_goes.sh
