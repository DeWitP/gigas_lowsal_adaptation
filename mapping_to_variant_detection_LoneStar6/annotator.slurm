#!/bin/bash
#SBATCH -J annotator
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p normal
#SBATCH -o annotator.o%j
#SBATCH -e annotator.e%j
#SBATCH -t 1:00:00
#SBATCH -A IBN21018
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pierre.de_wit@gu.se

java -Xmx32G -jar /scratch/05301/dewitp/dynamo/dedup_bams/GenomeAnalysisTK_3_6.jar -T VariantAnnotator \
    -V snps.vcf \
    -R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
    -A MappingQualityRankSumTest -A ReadPosRankSumTest \
    -o snps_annotated.vcf
    
