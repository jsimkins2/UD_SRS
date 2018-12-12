#!/bin/bash
# this script is designed to take raw VIIRS passes and use IPOPP and gpt software to process them to level 2 files

~/drl/tools/services.sh start
# first step is to remove all old files from IPOPP's ingest folder
rm *.PDS
rm *.h5

rm ~/drl/data/dsm/ingest/*
rm ~/drl/data/pub/gsfcdata/npp/viirs/level0/*
rm ~/drl/data/pub/gsfcdata/npp/viirs/level1/*
rm ~/viirs_ipopp_fix/*.h5*
rm ~/viirs_ipopp_fix/*.tif
rm ~/viirs_ipopp_fix/*.nc
rm ~/viirs_ipopp_fix/*.L2_1KM

# run the raw viirs file to create the PDS files
sh viirs_raw_L0.sh 
echo "finished the raw files, copying over to the ingest folder now"
# copy the files over to the ingest folder
scp *.PDS ~/drl/data/dsm/ingest/
scp *.h5 ~/drl/data/dsm/ingest/
scp *.dat ~/drl/data/dsm/ingest/
scp *.npp ~/drl/data/dsm/ingest/

# run the ingest_ipopp.sh
echo "running ipopp ingest now"
sh ~/drl/tools/ingest_ipopp.sh 

echo "ran ipopp, gonna take a nap for 10 before we check the files"
sleep 5m

#######################################################
# NEED TO RE-CODE THIS AT SOME POINT BUT NOT SURE OF A BETTER WAY JUST YET
# problem is that IPOPP will say it's finished before outputting files

# copy files into the viirs_ipopp_fix folder once they hit the level1 folder
# this if statement lists the files in the folder and when the word count list is greater than 30 we go forth with the code
# wc stands for word count, though in this fashion it's really used as a line counter from the list hence the -l argument
if [ `ls ~/drl/data/pub/gsfcdata/npp/viirs/level1/* | wc -l > 30` ]]
then
    scp ~/drl/data/pub/gsfcdata/npp/viirs/level1/* ~/viirs_ipopp_fix/
else
    echo "no files yet, sleeping for 5 minutes"
    sleep 5m
    if [ `ls ~/drl/data/pub/gsfcdata/npp/viirs/level1/* | wc -l > 30` ]]
    then
        scp ~/drl/data/pub/gsfcdata/npp/viirs/level1/* ~/viirs_ipopp_fix/
    else
        echo "no files yet, sleeping for 5 minutes"
        sleep 5m
        if [ `ls ~/drl/data/pub/gsfcdata/npp/viirs/level1/* | wc -l > 30` ]]
        then
            scp ~/drl/data/pub/gsfcdata/npp/viirs/level1/* ~/viirs_ipopp_fix/
        else
            echo "no files yet, sleeping for 5 minutes"
            sleep 5m
            if [ `ls ~/drl/data/pub/gsfcdata/npp/viirs/level1/* | wc -l > 30` ]]
            then
                scp ~/drl/data/pub/gsfcdata/npp/viirs/level1/* ~/viirs_ipopp_fix/
            else
                echo "no files yet, sleeping for 5 minutes"
                sleep 5m
                if [ `ls ~/drl/data/pub/gsfcdata/npp/viirs/level1/* | wc -l > 30` ]]
                then
                    scp ~/drl/data/pub/gsfcdata/npp/viirs/level1/* ~/viirs_ipopp_fix/
                else
                    echo "stilllll no files yet, sleeping one last time for 10 minutes before aborting mission"
                    sleep 10m
                    if [ `ls ~/drl/data/pub/gsfcdata/npp/viirs/level1/* | wc -l > 30` ]]
                    then
                        scp ~/drl/data/pub/gsfcdata/npp/viirs/level1/* ~/viirs_ipopp_fix/
                    fi
                fi
            fi
        fi
    fi
fi

# alright we now have the level1 files form IPOPP time to run l2gen
echo "got the level1 files, time to run viirs_l2.sh !"

sh ~/viirs_ipopp_fix/viirs_L2.sh 

echo "all done!"


~/drl/tools/services.sh stop




