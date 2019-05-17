#!/usr/bin/env bash

featureCounts -a $gtf -o $file $bam 

# 如果一个个样本单独计数，输出多个文件使用代码是：
for fn in {508..523}
do
featureCounts -T 5 -p -t exon -g gene_id  -a /public/reference/gtf/gencode/gencode.v25.annotation.gtf.gz -o $fn.counts.txt SRR1039$fn.bam
done

# 如果是批量样本的bam进行计数，使用代码是：
mkdir $wkd/align 
cd $wkd/align 
source activate rna

gtf="/public/reference/gtf/gencode/gencode.v25.annotation.gtf.gz"   
featureCounts -T 5 -p -t exon -g gene_id  -a $gtf -o  all.id.txt  *.bam  1>counts.id.log 2>&1 &

# https://blog.csdn.net/herokoking/article/details/78257714
# samtools view SRR1039516.subjunc.bam chr1:2350414-2352820|wc
less -SN all.id.txt


multiqc all.id.txt.summary
