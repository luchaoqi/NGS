# Sources:  
*  https://github.com/stephenturner/oneliners/blob/master/README.md#contents  
# Commands:
fastq sequences length distribution:  
```
zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'  
```
