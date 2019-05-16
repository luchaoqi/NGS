#!/usr/bin/env bash

# hisat2,subjunc,star,bwa,bowtie2 can be used to do alignment
# note need to create index at first


# subsample: cat for fq, zcat for gz
# if no need to keep .gz as file extension use basename $id ".gz"
ls /home/leon/rna-seq/project_Homo_sapiens/2_quality_control/trim_galore_output/*.fq| while read id; do ( cat $id|head -1000 > $(basename $id ));done

# bowtie2

bowtie2-build chr22.fa chr22
bowtie2 -p 10 -x chr22 -1 sample/pair.1.fq -2 sample/pair.2.fq -S sample.sam

# hisat2
hisat2-build ref.fa ref
hisat2 -p 10 -x ref -1 $fq1 -2 $fq2 -S temp.sam

# subjunc output bam file directly but for further processing purpose named as .sam file
subjunc -T 5 -i $index -r $fq1 -R $fq2 -o temp.sam

# bwa
bwa mem -t 5 -M $index $fq1 $fq2 > temp.sam


# alignment after trimming QC
ls *gz|cut -d"_" -f 1 |sort -u |while read id;do
ls -lh ${id}_1_val_1.fq.gz   ${id}_2_val_2.fq.gz 
hisat2 -p 10 -x /public/reference/index/hisat/hg38/genome -1 ${id}_1_val_1.fq.gz   -2 ${id}_2_val_2.fq.gz  -S ${id}.hisat.sam
subjunc -T 5  -i /public/reference/index/subread/hg38 -r ${id}_1_val_1.fq.gz -R ${id}_2_val_2.fq.gz -o ${id}.subjunc.sam
bowtie2 -p 10 -x /public/reference/index/bowtie/hg38  -1 ${id}_1_val_1.fq.gz   -2 ${id}_2_val_2.fq.gz  -S ${id}.bowtie.sam
bwa mem -t 5 -M  /public/reference/index/bwa/hg38   ${id}_1_val_1.fq.gz   ${id}_2_val_2.fq.gz > ${id}.bwa.sam
done


# sam2bam
ls *.sam|while read id; do ( samtools sort -O bam -@ 5 -o $(basename $id ".sam").bam $id );done

# index for bam
ls *.bam|xargs -i samtools index {}

samtools view temp.bam | less -S
samtools tview temp.bam
# can use igv instead of tview
samtools flagstat temp.bam



# QC comparison
ls *.bam |while read id ;do ( samtools flagstat -@ 1 $id >  $(basename ${id} ".bam").flagstat  );done

# report excel
# https://www.youtube.com/watch?v=Pmws2RmlJRU&list=PLKOVv6BeGBMQ2bmWE2DRD2NYk5jPh-huJ&index=12
# cat *.flagstat|awk '{print $1}'|paste - - - - - - - - - - - - - 

multiqc *.flagstat
