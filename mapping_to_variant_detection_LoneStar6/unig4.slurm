#!/bin/bash
#SBATCH -J gt
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p normal
#SBATCH -o gt.o%j
#SBATCH -e gt.e%j
#SBATCH -t 24:00:00
#SBATCH -A IBN21018
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pierre.de_wit@gu.se

java -Xmx64G -jar /scratch/05301/dewitp/dynamo/dedup_bams/GenomeAnalysisTK_3_6.jar -T UnifiedGenotyper \
-R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -nt 40 -nct 1 \
--genotype_likelihoods_model SNP \
-I F01.realigned.bam \
-I F02.realigned.bam \
-I F03.realigned.bam \
-I F04.realigned.bam \
-I F05.realigned.bam \
-I F06.realigned.bam \
-I F07.realigned.bam \
-I F08.realigned.bam \
-I F09.realigned.bam \
-I F10.realigned.bam \
-I F11.realigned.bam \
-I F12.realigned.bam \
-I F13.realigned.bam \
-I F14.realigned.bam \
-I F15.realigned.bam \
-I M01.realigned.bam \
-I M02.realigned.bam \
-I M03.realigned.bam \
-I M04.realigned.bam \
-I M05.realigned.bam \
-I M06.realigned.bam \
-I M07.realigned.bam \
-I M08.realigned.bam \
-I M09.realigned.bam \
-I M10.realigned.bam \
-I M11.realigned.bam \
-I M12.realigned.bam \
-I M13.realigned.bam \
-I M14.realigned.bam \
-I M15.realigned.bam \
-o parents_raw_SNPs.vcf