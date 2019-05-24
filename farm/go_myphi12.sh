#!/bin/bash

INPUT_FILE=$1
OUTPUT_FILE=$2
RUN=$3

echo " >>  $INPUT_FILE $OUTPUT_FILE $RUN "
cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/phi_analysis/pid/filter_clas12_v4.cxx .

#root -l<<EOF
#.L filter_clas12_v3.cxx++ 
#filter_clas12_v3("$INPUT_FILE","$OUTPUT_FILE",$RUN) 
#.q
#EOF


root -l -b -q "filter_clas12_v4.cxx+("\"$INPUT_FILE\"","\"$OUTPUT_FILE\"",$RUN)"
