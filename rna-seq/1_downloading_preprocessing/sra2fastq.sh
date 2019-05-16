# will dump reads in a number of "standard" fastq and fasta formats.
# just check fastq-dump --h
fastq-dump --gzip --split-3 -O ./ your sra path

# e.g.
ls /public/project/RNA/airway/sra/* | while read id; do ( nohup fastq-dump --gzip --split-3 -O ./ $id &); done

