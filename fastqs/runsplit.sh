#!/bin/bash

for dir in 1465-*; do
    cd $dir
    echo $dir
    sampleid=`cat sampleid.txt`
    for fq in `find . -name '*.fastq.gz'`; do
        mate=$(echo $fq | grep -oP 'WGSb_[0-9]+.rgFC' | grep -oP '[0-9]+')
        bsub -q short -W 12:00 -o $fq.split.log ../split_and_rename.py $sampleid $mate $fq
    done
    cd ..
done
