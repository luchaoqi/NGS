# Sources:  
*  https://github.com/stephenturner/oneliners/blob/master/README.md#contents  
*  https://mp.weixin.qq.com/s?__biz=MzU4NjU4ODQ2MQ==&mid=2247486039&idx=2&sn=887b456b6b762b75094eb70cfbad1ac6&chksm=fdf84215ca8fcb0308213948c32cf7a70449e2b5a0c2fcee4f6b80f650738ec1e00258934c96&mpshare=1&scene=1&srcid=&key=c8c9cb9453e0935031f4a12f7d3a5784c4026682a6cf81d6b04d42e00e93b2da51c4037221d45e851a0e91480113444778b6596784c8cbf3a64a28976ac19bc67a12a323c5af0283a5f16be6e2ec5038&ascene=1&uin=NzA3NTE3MTMz&devicetype=Windows+10&version=62060739&lang=en&pass_ticket=qsvmCxoFviqMD1c7jgwYAO3C7Uu8eH2USUAm%2F8ynoPIf5aNtlWgXD2Y1sNQLgw4I  
# Commands:
fastq sequences length distribution  => 得到fq文件中序列长度的分布:  
```
zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'  
```

reverse complement  => 反向互补:  

```
echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC'
```

fastq2fasta:  

```
zcat file.fastq.gz | paste - - - - | perl -ane 'print ">$F[0]\n$F[1]\n";' | gzip -c > file.fasta.gz
```

split a multifasta file into single ones with csplit => fasta按>拆分:  

```
csplit -z -q -n 4 -f test test.fa /\>/ {*}
```
```
awk '/^>/{s=++d".fa"} {print > s}' multi.fa
```

single line fasta to multi-line of 50 characters in each line => 单行fa变多行:  

```
awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/50;i++) {for (j=1;j<=50;j++) printf "%s", $(i*50 +j); print ""}}' file
```

multi-line fasta to one-line => 一个多行fa文件变单行:  

```
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' file.fa
```
```
cat file.fasta | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}'
```

bam2bed:  

```
samtools view file.bam | perl -F'\t' -ane '$strand=($F[1]&16)?"-":"+";$length=1;$tmp=$F[5];$tmp =~ s/(\d+)[MD]/$length+=$1/eg;print "$F[2]\t$F[3]\t".($F[3]+$length)."\t$F[0]\t0\t$strand\n";' > file.bed
```
bam2wig:  

```
samtools mpileup -BQ0 file.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > file.wig.gz
```
