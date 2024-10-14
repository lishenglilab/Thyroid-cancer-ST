#!/bin/bash

samplePath="/home/yvzeng/ST/"
samples="ATC-1 ATC-2 ATC-3 ATC-4 PTC-1 PTC-2 PTC-3 PTC-4 LPTC-1 LPTC-2 LPTC-3 LPTC-4 N-1 N-2 N-3 N-4"

for i in $samples
do
	spaceranger count --id=${i} \
                   --transcriptome=/opt/refdata/GRCh38-3.0.0 \
                   --fastqs=${samplePath}/${i} \
                   --sample=${i} \
                   --image=${samplePath}/${i}/${i}.tif \
                   --slide=V19J01-123 \
                   --area=A1
done
