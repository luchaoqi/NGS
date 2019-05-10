# Sources:  
*  https://github.com/stephenturner/oneliners/blob/master/README.md#contents  
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
