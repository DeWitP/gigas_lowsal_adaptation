#!/bin/bash

## Command log for analysis of the DynamO crossing experiment data, Oct-Nov 2023.
# Adults are located in $SCRATCH/dynamo/adults
# Embryos in $SCRATCH/dynamo/embryos

# Make sure that the right tools are installed:

mamba activate analysis

# mamba install angsd
# mamba install bowtie2
# mamba install cutadapt
# mamba install fastqc
# mamba install samtools
# mamba install -c conda-forge r-base=4.1.2
# mamba install -c bioconda trim-galore 
# mamba install -c bioconda bwa
# mamba install -c bioconda picard
# mamba install vcftools
# mamba install bcftools=1.10
# mamba install bedops
# mamba install jvarkit-wgscoverageplotter
# mamba install freebayes


# Indexing the genome sequence (From Penaloza et al 2021: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_902806645.1/)
export GENOME_FASTA=/scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna
export GENOME_DICT=/scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.dict 
bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA
bwa index $GENOME_FASTA

idev -A IBN21018
picard CreateSequenceDictionary R=$GENOME_FASTA O=$GENOME_DICT
exit

#Note - picard does not work on the login node, something with java.

###### ----------- ######
# STEP ONE. QC AND MAPPING 

# First, working with the adults. 

#Renaming the adult files from the machine names to their biological names:

./rename_adults.sh
# FastQC was already run by SciLife,
# Found very high quality, but some non-randomness in the first 9 bases, and also some Nextera adapter.
# So, removing 9 bases at 5' end and Nextera adapters with TrimGalore, keeping reads longer than 50 bases after trimming:

>trimclip
for file in *_R1.fastq.gz; do
echo "trim_galore -q 30 -a CTGTCTCTTATA -a2 CTGTCTCTTATA --paired --length 50 --clip_R1 9 --clip_R2 9 --retain_unpaired ${file/_R1.fastq.gz/}_R1.fastq.gz ${file/_R1.fastq.gz/}_R2.fastq.gz">>trimclip;
done

ls6_launcher_creator.py -j trimclip -n trimclip -t 1:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch trimclip.slurm

# Mapping the short read data to the genome:

export GENOME_FASTA=/scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna
            
for F in `ls *_R1_val_1.fq.gz`; 
	do base="$(basename $F _R1_val_1.fq.gz)"; 
	ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")
	echo "bwa mem -t 12 -R \"$ID\" $GENOME_FASTA "$base"_R1_val_1.fq.gz "$base"_R2_val_2.fq.gz >"$base".sam" >> BWAMEMMAP;
done
            
ls6_launcher_creator.py -j BWAMEMMAP -n BWAMEMMAP -t 4:00:00 -w 4 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch BWAMEMMAP.slurm

ls *sam | wc -l # should be the same number as number of fastq files

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
             ls *.sam >sams
             
# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:

>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

ls6_launcher_creator.py -j s2b -n s2b -t 2:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch s2b.slurm

ls *bam | wc -l  # should be the same number as number of trim files
rm *.sam
ls *.bam >bams   

###------------------------------------------------------
# Now, working with the embryos

cd ../embryos
./rename_embryos.sh

# FastQC data for embryos show same pattern as for adults. Good quality, but first 9 bases biased and some Nextera adapters.
# So, running the same TrimGalore command.

>trimclip
for file in *_R1.fastq.gz; do
echo "trim_galore -q 30 -a CTGTCTCTTATA -a2 CTGTCTCTTATA --paired --length 50 --clip_R1 9 --clip_R2 9 --retain_unpaired ${file/_R1.fastq.gz/}_R1.fastq.gz ${file/_R1.fastq.gz/}_R2.fastq.gz">>trimclip;
done

ls6_launcher_creator.py -j trimclip -n trimclip -t 1:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch trimclip.slurm

# Mapping the short read data to the genome:

export GENOME_FASTA=/scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna
            
for F in `ls *_R1_val_1.fq.gz`;
	do base="$(basename $F _R1_val_1.fq.gz)";
	ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")
	echo "bwa mem -t 12 -R \"$ID\" $GENOME_FASTA "$base"_R1_val_1.fq.gz "$base"_R2_val_2.fq.gz >"$base".sam" >> BWAMEMMAP;
done
            
ls6_launcher_creator.py -j BWAMEMMAP -n BWAMEMMAP -t 4:00:00 -w 4 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch BWAMEMMAP.slurm

ls *sam | wc -l # should be the same number as number of fastq files

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
             ls *.sam >sams
             
# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:

>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

ls6_launcher_creator.py -j s2b -n s2b -t 2:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch s2b.slurm

#### -------------------- #####
# STEP TWO. WORKING WITH THE BAM FILES

# Removing duplicate reads.
 
for F in `ls *.bam`;
	do base="$(basename $F .bam)"
	echo "picard MarkDuplicates INPUT="$base".bam OUTPUT="$base".dedup.bam METRICS_FILE="$base".metrics AS=true CREATE_INDEX=true" >> MD;
done
 
ls6_launcher_creator.py -j MD -n MD -t 1:00:00 -w 6 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch MD.slurm
             
# INDEL REALIGNMENT
             
sbatch targetcreator.slurm
             
for F in `ls *.dedup.bam`;
	do base="$(basename $F .dedup.bam)"
echo "java -Xmx32G -jar /scratch/05301/dewitp/dynamo/dedup_bams/GenomeAnalysisTK_3_6.jar -T IndelRealigner \
-R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
-targetIntervals realigner.intervals \
-I "$base".dedup.bam \
-o "$base".realigned.bam" >> realigner;
done 
             
ls6_launcher_creator.py -j realigner -n realigner -t 6:00:00 -w 4 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch realigner.slurm
 
 
## Merging the three embryo bams for each sample:

for F in `ls *_14.realigned.bam`;
	do base="$(basename $F _14.realigned.bam)"

	echo "samtools merge "$base"_14.merged.realigned.bam "$base"_14.realigned.bam "$base"_14_1.realigned.bam "$base"_14_2.realigned.bam" >> merger
done
for F in `ls *_24.realigned.bam`;
	do base="$(basename $F _24.realigned.bam)"

	echo "samtools merge "$base"_24.merged.realigned.bam "$base"_24.realigned.bam "$base"_24_1.realigned.bam "$base"_24_2.realigned.bam" >> merger
done

ls6_launcher_creator.py -j merger -n merger -t 6:00:00 -w 4 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch merger.slurm

for F in `ls *.merged.realigned.bam`;
	do
	echo "samtools index $F" >> index_merged;
done

ls6_launcher_creator.py -j index_merged -n index_merged -t 0:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch index_merged.slurm
  
  
# Investigating coverage of the merged embryo bams using wgscoverageplotter

for F in `ls *.merged.realigned.bam`;
	do base="$(basename $F .realigned.bam)"
echo "wgscoverageplotter.py --dimension 1500x500 -C -1 --clip -R ../ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna ./"$base".realigned.bam --min-contig-length 30mb --percentile median  > ./"$base"_cov.svg" >> covplot_merged;
done

ls6_launcher_creator.py -j covplot_merged -n covplot_merged -t 0:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch covplot_merged.slurm


# green line - mean, red line - median depth. Looks OK!

###----------------------------------###
# STEP THREE. GENOTYPING THE PARENTS

#Genotyping the parents with GATK:
             
sbatch unig4.slurm

### Hard filtering variants (according to https://gatk.broadinstitute.org/hc/en-us/articles/360035531112#2):

# First, selecting SNP-only variants:

java -Xmx32G -jar /scratch/05301/dewitp/dynamo/dedup_bams/GenomeAnalysisTK_3_6.jar -T SelectVariants \
    -V parents_raw_SNPs.vcf \
    -R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
    -selectType SNP \
    -o snps.vcf


java -Xmx32G -jar /scratch/05301/dewitp/dynamo/dedup_bams/GenomeAnalysisTK_3_6.jar -T VariantAnnotator \
    -V snps.vcf \
    -R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
    -A MappingQualityRankSumTest -A ReadPosRankSumTest \
    -o snps_annotated.vcf


# Then, applying the standard hard filters recommended by the Broad Institute:

java -Xmx32G -jar /scratch/05301/dewitp/dynamo/dedup_bams/GenomeAnalysisTK_3_6.jar -T VariantFiltration \
    -V snps_annotated.vcf \
    -R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
    -filter "QD < 2.0" --filterName "QD2" \
    -filter "QUAL < 30.0" --filterName "QUAL30" \
    -filter "SOR > 3.0" --filterName "SOR3" \
    -filter "FS > 60.0" --filterName "FS60" \
    -filter "MQ < 40.0" --filterName "MQ40" \
    -filter "MQRankSum < -12.5" --filterName "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum-8" \
    -o snps_filtered.vcf

#NOTE - Get warning for all loci where all individuals are non-ref homozygotes (Missing MQRankSum and ReadPosRankSum) - will filter all of these loci out below.

#Finally, save all the SNPS that have passed the VQSR filter into a new vcf file:

grep "PASS\|^#" snps_filtered.vcf > PASSING_SNPS.vcf

# Let's do some vcftools filtering too! The below keeps only biallelic loci, no missing data.

vcftools --vcf PASSING_SNPS.vcf --remove-filtered-all --max-missing 1  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out biallelic_nomissing2

# After filtering, kept 24838954 out of a possible 45513499 Sites

# Turns out that for some of these, all individuals are non-ref homozygotes. Adding the --max-non-ref-ac flag to remove these:

vcftools --vcf PASSING_SNPS.vcf --remove-filtered-all --max-missing 1  --min-alleles 2 --max-alleles 2 --max-non-ref-ac 29 --recode --recode-INFO-all --out biallelic_nomissing2

# After filtering, kept 22132741 out of a possible 45513499 Sites


#Finally, output all kinds of statistics with vcftools:
vcftools --vcf biallelic_nomissing2.recode.vcf --out biallelic_nomissing2 --het
vcftools --vcf biallelic_nomissing2.recode.vcf --out biallelic_nomissing2 --window-pi 10000
vcftools --vcf biallelic_nomissing2.recode.vcf --out biallelic_nomissing2 --relatedness

# Note about "relatedness" from the vcftools manual:
# This option is used to calculate and output a relatedness statistic based on the method of Yang et al, Nature Genetics 2010 
# (doi:10.1038/ng.608). Specifically, calculate the unadjusted Ajk statistic. Expectation of Ajk is zero for individuals 
# within a populations, and one for an individual with themselves. The output file has the suffix ".relatedness".

vcftools --vcf biallelic_nomissing2.recode.vcf --out biallelic_nomissing2 --TajimaD 10000

#Test with long runs of homozygosity for one chromosome, indicator of recent bottleneck:
vcftools --vcf biallelic_nomissing2.recode.vcf --out biallelic_nomissing2 --LROH --chr LR761634.1
vcftools --vcf biallelic_nomissing2.recode.vcf --out biallelic_nomissing2_chr6 --LROH --chr LR761639.1
vcftools --vcf biallelic_nomissing2.recode.vcf --out biallelic_nomissing2_chr9 --LROH --chr LR761642.1

###-------------------------####
#STEP FOUR. EXTRACTING ALLELE READ COUNTS FROM THE EMBRYOS

 		# Extracting the allele counts from the juvenile data using samtools mpileup:
		# Syntax: samtools mpileup [-EB] [-C capQcoef] [-r reg] [-f in.fa] [-l list] [-Q minBaseQ] [-q minMapQ] in.bam [in2.bam [...]] 
		
		# First, extracting the list of SNP loci as a bed file from the vcf file
		vcf2bed < biallelic_nomissing2.recode.vcf > SNPs.bed
		cut -f 1-3 SNPs.bed > locilist.bed

		
		echo "samtools mpileup -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -l locilist.bed --ignore-RG -Q 30 -q 30 -o embryo_read_counts.mpileup -b mergedbamlist.txt" > mpileup
		
		ls6_launcher_creator.py -j mpileup -n mpileup -t 4:00:00 -w 1 -a IBN21018 -q normal -e pierre.de_wit@gu.se
		sbatch mpileup.slurm
		
		# samtools mpileup takes too long! Trying freebayes instead:

# Zip and index the vcf file for freebayes use as reference:
bgzip biallelic_nomissing2.recode.vcf
tabix biallelic_nomissing2.recode.vcf.gz

# First, have to homogenize the read groups in the merged bam files, otherwise FreeBayes separates them in the vcf file.

for F in `ls *.merged.realigned.bam`;
	do base="$(basename $F .merged.realigned.bam)"
echo "picard AddOrReplaceReadGroups INPUT=$base.merged.realigned.bam OUTPUT=$base.RGmerged.bam RGID=$base RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$base" >>RGmerge;
done

ls6_launcher_creator.py -j RGmerge -n RGmerge -t 1:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch RGmerge.slurm

#Indexing the new bam files:

for F in `ls *.RGmerged.bam`;
	do
echo "samtools index $F" >> index_RGmerged;
done

ls6_launcher_creator.py -j index_RGmerged -n index_RGmerged -t 0:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch index_RGmerged.slurm

#Running freebayes in pooled-seq data mode with unknown ploidy, using the parent vcf to extract positions (instructions here:https://github.com/freebayes/freebayes):

# The following commented out lines are a trial and error test of the freebayes commands. Keeping it here just for posteriority.

# # for F in `ls *.RGmerged.bam`;
# # 	do base="$(basename $F .RGmerged.bam)"
# # echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz $F >$base.vcf" >>fb;
# # done
# # 
# # ls6_launcher_creator.py -j fb -n fb -t 12:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# # sbatch fb.slurm
# 
# # The above timed out after 12 hours, pulled out a lot more loci than in the parent vcf. Maybe should add the --only-use-input-alleles flag.
# 
# # echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --only-use-input-alleles -L bamlist.txt >embryos_combined.vcf" > fb2
# ls6_launcher_creator.py -j fb2 -n fb2 -t 24:00:00 -w 1 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# sbatch fb2.slurm
# # Above only took 20 mins to run, produced an empty vcf. 
# 
# 
# # Trying the locilist instead, with the --targets flag, see if that works better...
# # Tried a joint calling with all embryo files, didn't work, so instead doing a separate calling for each embryo samples, and then merging the vcf files.
# 
# for F in `ls *.RGmerged.bam`;
# 	do base="$(basename $F .RGmerged.bam)"
# echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous --targets locilist.bed $F >$base.vcf" >> fb3;
# done
# ls6_launcher_creator.py -j fb3 -n fb3 -t 12:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# sbatch fb3.slurm
# 
# #Output looks good, but could only do about 6.5 M SNPs in 12 hours. Testing to divide the locilist into 4 parts.
# 
# for F in `ls *.RGmerged.bam`;
# 	do base="$(basename $F .RGmerged.bam)"
# echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous --targets locilist1.bed $F >$base.1.vcf" >> fb4_1;
# echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous --targets locilist2.bed $F >$base.2.vcf" >> fb4_2;
# echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous --targets locilist3.bed $F >$base.3.vcf" >> fb4_3;
# echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous --targets locilist4.bed $F >$base.4.vcf" >> fb4_4;
# done
# ls6_launcher_creator.py -j fb4_1 -n fb4_1 -t 12:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# ls6_launcher_creator.py -j fb4_2 -n fb4_2 -t 12:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# ls6_launcher_creator.py -j fb4_3 -n fb4_3 -t 12:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# ls6_launcher_creator.py -j fb4_4 -n fb4_4 -t 12:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# sbatch fb4_1.slurm
# sbatch fb4_2.slurm
# sbatch fb4_3.slurm
# sbatch fb4_4.slurm
# 
# #vcf files contain ca. 1.8 M SNPs, should be 6. Suppose that uncovered and monomorphic ones don't get output. Trying again with the adult vcf, but only one sample 
# 
# echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --only-use-input-alleles F9_M9_24.RGmerged.bam >F9_M9_24.vcf" > fbtest
# ls6_launcher_creator.py -j fbtest -n fbtest -t 8:00:00 -w 1 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# sbatch fbtest.slurm
# 
# # Doesn't work, outputs empty file. How about adding the locilist as well as the input vcf??
# 
# echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --targets locilist.bed --only-use-input-alleles -L bamlist.txt >embryos_combined.vcf" > fb5
# ls6_launcher_creator.py -j fb5 -n fb5 -t 12:00:00 -w 1 -a IBN21018 -q normal -e pierre.de_wit@gu.se
# sbatch fb5.slurm
# 
# ## Looked like it was working in interactive mode!!!  But only output 200 K SNPs in 12 hours...
# Trying again the sample-specific calling, now with the input vcf and also giving it 24 hours to run.

for F in `ls *.RGmerged.bam`;
	do base="$(basename $F .RGmerged.bam)"
echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --targets locilist_top10M.bed --only-use-input-alleles $F >$base.vcf" >> fb7_top;
done
ls6_launcher_creator.py -j fb7_top -n fb7_top -t 24:00:00 -w 4 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch fb7_top.slurm

for F in `ls *.RGmerged.bam`;
	do base="$(basename $F .RGmerged.bam)"
echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --targets locilist_bottom12M.bed --only-use-input-alleles $F >$base.2.vcf" >> fb7_bottom;
done
ls6_launcher_creator.py -j fb7_bottom -n fb7_bottom -t 24:00:00 -w 4 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch fb7_bottom.slurm

# The above works for 96 of the commands, but the other 48 produces an error "signal11" which aborts freebayes after 1-2 hours. All 48 are in the top 10M loci. Why???
# Tried running with w=1, in case RAM usage was the issue, but got same fault.
# Apparently, this problem has been described by other people too in freebayes 1.3.6: https://github.com/freebayes/freebayes/issues/390
# There is no newer version that works with conda.

## NOTE: All those files die on the same place: CADCXH010000226.1	54748
# That locus is the last of the satellite contigs in the assembly! AFter that starts the real chromosomes. What is going on?
# Creating a new locilist_top10M which excludes all the satellite contigs, try to run using that instead.

freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --only-use-input-alleles --targets locilist_top10M_nosatellites.bed F7_M7_14.RGmerged.bam > F7_M7_14.vcf

# WORKED! Took 10 h to run. So, let's just skip the satellite contigs and focus on the main chromosomes.

for F in `ls *.RGmerged.bam`;
	do base="$(basename $F .RGmerged.bam)"
echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --targets locilist_top10M_nosatellites.bed --only-use-input-alleles $F >$base.vcf" >> fb8_top;
done
ls6_launcher_creator.py -j fb8_top -n fb8_top -t 24:00:00 -w 4 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch fb8_top.slurm

# Then, do the satellites separately:
for F in `ls *.RGmerged.bam`;
	do base="$(basename $F .RGmerged.bam)"
echo "freebayes -f /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna -F 0.01 -C 1 --pooled-continuous -@ ../vcf/biallelic_nomissing2.recode.vcf.gz --targets locilist_satellites.bed --only-use-input-alleles $F >$base.3.vcf" >> fb9_sat;
done
ls6_launcher_creator.py -j fb9_sat -n fb9_sat -t 4:00:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch fb9_sat.slurm

# Merging all the bottom 12M files into a multi-sample vcf:
for F in `ls *.2.vcf`;
	do echo "bgzip $F " >> bgzip_vcfs_2
echo "tabix $F.gz " >> tabix_vcfs_2;
done
ls6_launcher_creator.py -j bgzip_vcfs_2 -n bgzip_vcfs_2 -t 1:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
ls6_launcher_creator.py -j tabix_vcfs_2 -n tabix_vcfs_2 -t 0:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch bgzip_vcfs_2.slurm
sbatch tabix_vcfs_2.slurm

echo "#!/bin/bash
#SBATCH -J merge_2
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p normal
#SBATCH -o merge_2.o%j
#SBATCH -e merge_2.e%j
#SBATCH -t 4:00:00
#SBATCH -A IBN21018
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pierre.de_wit@gu.se

bcftools merge F10_M10_14.2.vcf.gz F12_M10_14.2.vcf.gz F14_M14_14.2.vcf.gz F1_M2_24.2.vcf.gz F3_M3_14.2.vcf.gz F8_M7_14.2.vcf.gz F10_M10_24.2.vcf.gz F12_M10_24.2.vcf.gz F14_M14_24.2.vcf.gz F1_M3_24.2.vcf.gz F6_M4_14.2.vcf.gz F8_M7_24.2.vcf.gz F10_M11_14.2.vcf.gz F12_M11_14.2.vcf.gz F14_M15_14.2.vcf.gz F2_M1_14.2.vcf.gz F6_M4_24.2.vcf.gz F8_M8_14.2.vcf.gz F10_M11_24.2.vcf.gz F12_M11_24.2.vcf.gz F14_M15_24.2.vcf.gz F2_M1_24.2.vcf.gz F6_M5_24.2.vcf.gz F8_M8_24.2.vcf.gz F10_M12_14.2.vcf.gz F12_M12_14.2.vcf.gz F15_M13_14.2.vcf.gz F2_M2_14.2.vcf.gz F6_M6_14.2.vcf.gz F8_M9_14.2.vcf.gz F10_M12_24.2.vcf.gz F12_M12_24.2.vcf.gz F15_M13_24.2.vcf.gz F2_M2_24.2.vcf.gz F6_M6_24.2.vcf.gz F8_M9_24.2.vcf.gz F11_M10_14.2.vcf.gz F13_M13_14.2.vcf.gz F15_M14_14.2.vcf.gz F2_M3_14.2.vcf.gz F7_M7_14.2.vcf.gz F9_M7_14.2.vcf.gz F11_M10_24.2.vcf.gz F13_M13_24.2.vcf.gz F15_M14_24.2.vcf.gz F2_M3_24.2.vcf.gz F7_M7_24.2.vcf.gz F9_M7_24.2.vcf.gz F11_M11_14.2.vcf.gz F13_M14_14.2.vcf.gz F15_M15_14.2.vcf.gz F3_M1_14.2.vcf.gz F7_M8_14.2.vcf.gz F9_M8_14.2.vcf.gz F11_M11_24.2.vcf.gz F13_M14_24.2.vcf.gz F15_M15_24.2.vcf.gz F3_M1_24.2.vcf.gz F7_M8_24.2.vcf.gz F9_M8_24.2.vcf.gz F11_M12_14.2.vcf.gz F13_M15_24.2.vcf.gz F1_M1_14.2.vcf.gz F3_M2_14.2.vcf.gz F7_M9_14.2.vcf.gz F9_M9_14.2.vcf.gz F11_M12_24.2.vcf.gz F14_M13_14.2.vcf.gz F1_M1_24.2.vcf.gz F3_M2_24.2.vcf.gz F7_M9_24.2.vcf.gz F9_M9_24.2.vcf.gz -o last12Mloci.vcf.gz -O z" > merge_2.slurm

sbatch merge_2.slurm
tabix last12Mloci.vcf.gz

# Merging all the top 10M files into a multi-sample vcf:
for F in `ls *4.vcf`;
	do echo "bgzip $F " >> bgzip_vcfs_1
echo "tabix $F.gz " >> tabix_vcfs_1;
done
ls6_launcher_creator.py -j bgzip_vcfs_1 -n bgzip_vcfs_1 -t 1:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
ls6_launcher_creator.py -j tabix_vcfs_1 -n tabix_vcfs_1 -t 0:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch bgzip_vcfs_1.slurm
sbatch tabix_vcfs_1.slurm

echo "#!/bin/bash
#SBATCH -J merge_1
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p normal
#SBATCH -o merge_1.o%j
#SBATCH -e merge_1.e%j
#SBATCH -t 4:00:00
#SBATCH -A IBN21018
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pierre.de_wit@gu.se

bcftools merge F10_M10_14.vcf.gz F12_M10_14.vcf.gz F14_M14_14.vcf.gz F1_M2_24.vcf.gz F3_M3_14.vcf.gz F8_M7_14.vcf.gz F10_M10_24.vcf.gz F12_M10_24.vcf.gz F14_M14_24.vcf.gz F1_M3_24.vcf.gz F6_M4_14.vcf.gz F8_M7_24.vcf.gz F10_M11_14.vcf.gz F12_M11_14.vcf.gz F14_M15_14.vcf.gz F2_M1_14.vcf.gz F6_M4_24.vcf.gz F8_M8_14.vcf.gz F10_M11_24.vcf.gz F12_M11_24.vcf.gz F14_M15_24.vcf.gz F2_M1_24.vcf.gz F6_M5_24.vcf.gz F8_M8_24.vcf.gz F10_M12_14.vcf.gz F12_M12_14.vcf.gz F15_M13_14.vcf.gz F2_M2_14.vcf.gz F6_M6_14.vcf.gz F8_M9_14.vcf.gz F10_M12_24.vcf.gz F12_M12_24.vcf.gz F15_M13_24.vcf.gz F2_M2_24.vcf.gz F6_M6_24.vcf.gz F8_M9_24.vcf.gz F11_M10_14.vcf.gz F13_M13_14.vcf.gz F15_M14_14.vcf.gz F2_M3_14.vcf.gz F7_M7_14.vcf.gz F9_M7_14.vcf.gz F11_M10_24.vcf.gz F13_M13_24.vcf.gz F15_M14_24.vcf.gz F2_M3_24.vcf.gz F7_M7_24.vcf.gz F9_M7_24.vcf.gz F11_M11_14.vcf.gz F13_M14_14.vcf.gz F15_M15_14.vcf.gz F3_M1_14.vcf.gz F7_M8_14.vcf.gz F9_M8_14.vcf.gz F11_M11_24.vcf.gz F13_M14_24.vcf.gz F15_M15_24.vcf.gz F3_M1_24.vcf.gz F7_M8_24.vcf.gz F9_M8_24.vcf.gz F11_M12_14.vcf.gz F13_M15_24.vcf.gz F1_M1_14.vcf.gz F3_M2_14.vcf.gz F7_M9_14.vcf.gz F9_M9_14.vcf.gz F11_M12_24.vcf.gz F14_M13_14.vcf.gz F1_M1_24.vcf.gz F3_M2_24.vcf.gz F7_M9_24.vcf.gz F9_M9_24.vcf.gz -o top10Mloci.vcf.gz -O z" > merge_1.slurm
sbatch merge_1.slurm
tabix top10Mloci.vcf.gz

# Merging all the satellite files into a multi-sample vcf:
for F in `ls *.3.vcf`;
	do echo "bgzip $F " >> bgzip_vcfs_3
echo "tabix $F.gz " >> tabix_vcfs_3;
done
ls6_launcher_creator.py -j bgzip_vcfs_3 -n bgzip_vcfs_3 -t 1:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
ls6_launcher_creator.py -j tabix_vcfs_3 -n tabix_vcfs_3 -t 0:30:00 -w 12 -a IBN21018 -q normal -e pierre.de_wit@gu.se
sbatch bgzip_vcfs_3.slurm
sbatch tabix_vcfs_3.slurm

echo "#!/bin/bash
#SBATCH -J merge_3
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p normal
#SBATCH -o merge_3.o%j
#SBATCH -e merge_3.e%j
#SBATCH -t 4:00:00
#SBATCH -A IBN21018
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pierre.de_wit@gu.se

bcftools merge F10_M10_14.3.vcf.gz F12_M10_14.3.vcf.gz F14_M14_14.3.vcf.gz F1_M2_24.3.vcf.gz F3_M3_14.3.vcf.gz F8_M7_14.3.vcf.gz F10_M10_24.3.vcf.gz F12_M10_24.3.vcf.gz F14_M14_24.3.vcf.gz F1_M3_24.3.vcf.gz F6_M4_14.3.vcf.gz F8_M7_24.3.vcf.gz F10_M11_14.3.vcf.gz F12_M11_14.3.vcf.gz F14_M15_14.3.vcf.gz F2_M1_14.3.vcf.gz F6_M4_24.3.vcf.gz F8_M8_14.3.vcf.gz F10_M11_24.3.vcf.gz F12_M11_24.3.vcf.gz F14_M15_24.3.vcf.gz F2_M1_24.3.vcf.gz F6_M5_24.3.vcf.gz F8_M8_24.3.vcf.gz F10_M12_14.3.vcf.gz F12_M12_14.3.vcf.gz F15_M13_14.3.vcf.gz F2_M2_14.3.vcf.gz F6_M6_14.3.vcf.gz F8_M9_14.3.vcf.gz F10_M12_24.3.vcf.gz F12_M12_24.3.vcf.gz F15_M13_24.3.vcf.gz F2_M2_24.3.vcf.gz F6_M6_24.3.vcf.gz F8_M9_24.3.vcf.gz F11_M10_14.3.vcf.gz F13_M13_14.3.vcf.gz F15_M14_14.3.vcf.gz F2_M3_14.3.vcf.gz F7_M7_14.3.vcf.gz F9_M7_14.3.vcf.gz F11_M10_24.3.vcf.gz F13_M13_24.3.vcf.gz F15_M14_24.3.vcf.gz F2_M3_24.3.vcf.gz F7_M7_24.3.vcf.gz F9_M7_24.3.vcf.gz F11_M11_14.3.vcf.gz F13_M14_14.3.vcf.gz F15_M15_14.3.vcf.gz F3_M1_14.3.vcf.gz F7_M8_14.3.vcf.gz F9_M8_14.3.vcf.gz F11_M11_24.3.vcf.gz F13_M14_24.3.vcf.gz F15_M15_24.3.vcf.gz F3_M1_24.3.vcf.gz F7_M8_24.3.vcf.gz F9_M8_24.3.vcf.gz F11_M12_14.3.vcf.gz F13_M15_24.3.vcf.gz F1_M1_14.3.vcf.gz F3_M2_14.3.vcf.gz F7_M9_14.3.vcf.gz F9_M9_14.3.vcf.gz F11_M12_24.3.vcf.gz F14_M13_14.3.vcf.gz F1_M1_24.3.vcf.gz F3_M2_24.3.vcf.gz F7_M9_24.3.vcf.gz F9_M9_24.3.vcf.gz -o satellites.vcf.gz -O z" > merge_3.slurm
sbatch merge_3.slurm
tabix satellites.vcf.gz

# FInally, concatenating the three merged files into one:

bcftools concat top10Mloci.vcf.gz last12Mloci.vcf.gz satellites.vcf.gz -o all_embryos_merged.vcf.gz -O z
tabix all_embryos_merged.vcf.gz

#Extracting allele depth from the embryo file using GATK:

# Testing with a small subset, 100 KB at the start of one chrom:
bcftools view locilistvcfs/last12Mloci.vcf.gz -r LR761639.1:1-100000 > first100KBtest/embryos.vcf

java -jar GenomeAnalysisTK_3_6.jar -T VariantsToTable \
-R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
-V embryos.vcf \
-F CHROM -F POS -F REF -F ALT -GF AD \
-o embryos.txt 

# And genotypes from the parent file:

java -jar GenomeAnalysisTK_3_6.jar -T VariantsToTable \
-R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
-V parents.vcf \
-F CHROM -F POS -F REF -F ALT -GF GT \
-o parents_genos.txt 

# Or maybe get the alternative allele counts instead:

vcftools --vcf parents.vcf --012

# Now, extracting data from the full dataset:

java -jar GenomeAnalysisTK_3_6.jar -T VariantsToTable \
-R /scratch/05301/dewitp/dynamo/ref/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
-V all_embryos_merged.vcf.gz \
-F CHROM -F POS -F REF -F ALT -GF AD \
-o embryos.txt 

vcftools --gzvcf biallelic_nomissing2.recode.vcf.gz --012 --out parents

# Finally, modify the 012 output a bit, to also get the individual and loci names into the same file as the genotype data:

cut -f 1 --complement parents.012 > testparents.txt
paste parents.012.indv testparents.txt > testparents2.txt
./transpose.pl #Input and output file names are hardcoded in the script

# Then, open the parents.012.pos file in nano and add a header line: CHROM	POS
# Re-save as parents_chrompos.txt

paste parents_chrompos.txt testparents2_transposed.txt > parents.txt

rm testparents*
rm parents_chrompos.txt

# Copy parents.txt and embryos.txt to another computer for analyses in R!!

 