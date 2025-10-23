#!/bin/bash
# OptiType batch processing script
# Loops through paired-end FASTQ files and runs OptiType on each sample

for R1 in SRA_Input/*_1.fastq.gz; do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    echo "Processing sample: $SAMPLE"
    OptiTypePipeline.py -i "$R1" "$R2" -r -o "SRA_Output/$SAMPLE/"
done

echo "All samples processed successfully."
