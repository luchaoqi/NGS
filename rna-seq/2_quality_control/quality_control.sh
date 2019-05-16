#!/usr/bin/env bash
# batch processing all files *.gz in ~/ncbi/public/sra/

# way to learn xargs
ls ~/ncbi/public/sra/*gz | xargs fastqc -t 10

# but usually just use fastqc directly
fastqc *.fq

# overview

multiqc ./


outputdir='/home/leon/rna-seq/project_Homo_sapiens/2_quality_control/trim_galore_output'
fq1='sample_data/frag180.1.fq'
fq2='sample_data/frag180.2.fq'

# patch processing all files, of note, just feed in fastq.gz file is also fine

ls /home/leon/rna-seq/project_Homo_sapiens/2_quality_control/sample_data/*1.fq > 1
ls /home/leon/rna-seq/project_Homo_sapiens/2_quality_control/sample_data/*2.fq > 2
paste 1 2 > config
cat config | while read id
#cat $1|while read id : bash quality_control.sh config
#cat $config_file |while read id
do
	arr=($id)
	fq1=${arr[0]}
	fq2=${arr[1]}
	nohup	trim_galore -q 25 --phred33 --length 30 -e 0.1 --stringency 3 --paired -o $outputdir $fq1 $fq2	&
done

cd $outputdir
conda activate python2
multiqc ./
