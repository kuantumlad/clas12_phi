#!/bin/bash

data_dir=/volatile/clas12/tylern/root/cooked_5.7.4/hipo2root/
sub_dir=/w/hallb-scifs17exp/clas12/bclary/CLAS12/phi_analysis/farm/
run=4013
fnum=0

total_files=1000 #`ls -l ${data_dir}/*${run}*.root | wc -l`
analysis_status="analysis"
priority_status="debug"
status=$analysis_status

echo " Number of runs to process for run 4013 " ${run}

while [ $fnum -lt $total_files ]
do
    echo " Creating jsub file job " ${fnum}
    
    if [ $fnum -gt 1000 ]
    then
	status=$analysis_status
    fi
    echo ${status}

    echo "JOBNAME: clas12PhiClary"${fnum} >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "OS: centos7" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub
    echo "TRACK: ${status}" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "PROJECT: clas12" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "MEMORY: 8096 MB" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "INPUT_FILES: ${data_dir}out_clas_00${run}.evio.${fnum}.root" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "OTHER_FILES: /w/hallb-scifs17exp/clas12/bclary/CLAS12/phi_analysis/farm/go_myphi12.sh" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "COMMAND: ./go_myphi12.sh out_clas_00${run}.evio.${fnum}.root filter_out_clas12_00${run}.${fnum}.root ${run}" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "OUTPUT_DATA: filter_out_clas12_00${run}.${fnum}.root" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 
    echo "OUTPUT_TEMPLATE: /lustre/expphy/volatile/clas12/bclary/clas12RunData/r00${run}/filter_out_clas12_00${run}.${fnum}.root" >> ${sub_dir}sub_data_phi_r${run}_job${fnum}.jsub 

    fnum=$((fnum+1))
done


i=0
while [ $i -lt $total_files ]
do
    echo " Submitting file ${i} to farm "
    jsub ${sub_dir}sub_data_phi_r${run}_job${i}.jsub 
    i=$((i+1))
done
