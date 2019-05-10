#!/usr/bin/env bash

# hisat2,subjunc,star,bwa,bowtie2 can be used to do alignment

# note need to create index at first
# subsample: cat for fq, zcat for gz
# if no need to keep .gz as file extension use basename $id ".gz"
ls /home/leon/rna-seq/project_Homo_sapiens/2_quality_control/trim_galore_output/*.fq| while read id; do ( cat $id|head -1000 > $(basename $id ));done

# bowtie2

bowtie2-build chr22.fa chr22
bowtie2 -x chr22 -1 sample/pair.1.fq -2 sample/pair.2.fq -S sample.sam


# hisat2
hisat2-build ref.fa ref
hisat2 -p 10 -x ref -1 $fq1 -2 $fq2 -S temp.sam

# subjunc
# build index at first!
subjunc -T 5 -i $index -r $fq1 -R $fq2 -o temp.sam

# bwa
bwa mem -t 5 -M $index $fq1 $fq2 > temp.sam


