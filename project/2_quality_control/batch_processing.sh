#!/usr/bin/env bash
# batch processing all files *.gz in ~/ncbi/public/sra/

# way to learn xargs
ls ~/ncbi/public/sra/*gz | xargs fastqc -t 10

# but usually just use fastqc directly
fastqc *.fq

# overview

multiqc ./
