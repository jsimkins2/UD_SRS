
#!/bin/sh
HOST='ftp.avl.class.noaa.gov'
USER='anonymous'
PASSWD='password'




for FILENAME in `cat listing001.txt | awk '{print $9}'`; do
ftp -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASSWD
binary
prompt no
cd /184029/7994388239/001
get $FILENAME
quit
END_SCRIPT
echo $FILENAME
done



