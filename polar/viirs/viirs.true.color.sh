# get ancillary data
getanc.py $1 > anc.txt

# remove 1st line of ancillary data which just says VIIRS
sed '1d' anc.txt > tmpfile; mv tmpfile anc.txt

# remove the file that getanc.py automatically makes with that ugly name
rm "$1.anc"
# grab the current directory for R
cwd=$(pwd)

# Use seadas l1brsgen to make the true color image and save as hdf file
l1brsgen ifile=$1 geofile=$2 ofile="$3.hdf" subsamp=8 atmocor=1 sline=3 resolution=1000 oformat=1 oformat_depth=24bit par=anc.txt

# run R and pass along the arguments necessary for jpeg output
Rscript hdf2jpeg.R ofile=$3 path=$cwd

# rotate the jpeg 180 degrees since it's always upside down for some reason
sips -r 180 "$3.jpeg"

# open up the true color jpeg!
open "$3.jpeg"
