#!/bin/bash

samplePath="/home/yvzeng/ST/SC"
samples="PTC-2 PTC-3 PTC-4 LPTC-2 LPTC-3 LPTC-4 N-2 N-3 N-4"

for i in $samples
do
	cellranger count --id=${i} --fastqs=${samplePath}/${i} --sample=${i} --transcriptome=/opt/ref/refdata-gex-GRCh38-2020-A
done
