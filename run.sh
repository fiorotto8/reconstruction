#!/bin/bash
#ped_run=$(awk -F ',' '{ if ($2 ~ /PED/ || $21==1) {print $1}}' RunlogfilteredbyDescription-data-2023-11-07_17_16_55.csv)
#ped_run=$(awk -F ',' '{ if ($2 ~ /PED/ || $21==1) {print $1}}' allRun.csv)
ped_run=$(awk -F ',' '{ if ($2 ~ /PED/ || $21==1) {print $1}}' ped_MANGO_NID.csv)

#original piece all the files
count=0
for run in $ped_run; do
    # Perform actions on each number, for example, print it
    echo "Processing number: $run"
    python3 reconstruction.py configFile_MANGOtemp.txt -r ${run} -t ../PedLIME --pdir ./output_plots -d ./output_files;
    count=$((1+$count))
done
echo $count
