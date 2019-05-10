#!/usr/bin/env bash

#dummy way to download
#for ((i=677;i<=680;i++)) ;do wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP029/SRP029245/SRR957$i/SRR957$i.sra;done
#sratools to download
cat SRR_Acc_List.txt | while read id; do ( prefetch $id & );done
# to kill the process
ps -ef | grep prefetch | awk '{print $2}'| while read id; do kill $id; done
