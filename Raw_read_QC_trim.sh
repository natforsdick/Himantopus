#!/bin/sh

#  Raw_read_QC_trim.sh
#  Created by Natalie Forsdick on 26/07/17.
#
# Here we do the first process following receipt of whole-genome short-read sequence data for kakī and pied stilts - 
# quality checking, adapter removal, trimming and deduplication.
# Programs: Trimmomatic-0.35, FastQC, condetri-v2.3.

# Where are the programs
fastqc=/usr/local/FastQC/fastqc  # new FASTQC as at 5/9/19 /usr/local/fastqc_v0.11.5/FastQC/fastqc
trimmomatic=/usr/local/Trimmomatic-0.35/trimmomatic-pe
condetri=/usr/local/condetri-v2.3/condetri.pl
dedup=/usr/local/condetri-v2.3/filterPCRdupl.pl

# Where is the data/output directory
datadir=/home/natalie/All_Raw/
home=/scratch/natalie/All_Trimmed/
now=$(date)

# Give it the samplist
samplist="LIST_OF_SAMP_NAMES"

# Tell it what the pairs are
fq1=_fq1.fastq
fq2=_fq2.fastq

##########
# 1) Unzip fastq files

if [ -e $datadir${samp}${fq1}.gz ] || [ -e $datadir${samp}${fq2}.gz ]; then
echo "Unzipping fastq files at $(date)"
gunzip -c $datadir${samp}${fq1}.gz > $home/${samp}${fq1}
gunzip -c $datadir${samp}${fq2}.gz > $home/${samp}${fq2}
fi

##########
# 2) Raw quality check

for samp in $samplist
do
echo "Checking raw quality of ${samp} at ${now}"
$fastqc ${samp}${fq1} 
$fastqc ${samp}${fq2}

done

##########
# 3) Adapter trimming only in Trimmomatic using the adapter reference file TruSeq3-PE-2.fa, and following this, remove reads < 50 bp.
# Output files will have the suffixes: _adaptertrimmed*.fastq, _adapterunpaired*.fastq.

for samp in $samplist
do
echo "Trimming adapters for ${samp} at ${now}"
$trimmomatic -threads 30 -phred33 $datadir${samp}${fq1} $datadir${samp}${fq2} ${samp}_adaptertrimmed1.fastq ${samp}_adapterunpaired1.fastq ${samp}_adaptertrimmed2.fastq ${samp}_adapterunpaired2.fastq ILLUMINACLIP:/usr/local/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:1:30:10 MINLEN:50
## updated Trimmomatic as at 5/9/19 - /usr/local/Trimmomatic-0.38/trimmomatic-pe

##########
# Optional FastQC between trimming steps.

#echo "Checking adapter trimmed quality of ${samp} at ${now}"
#$fastqc ${samp}_adaptertrimmed1.fastq 
#$fastqc ${samp}_adaptertrimmed2.fastq

done


#########################
# 4) Quality trimming reads with condetri.
# Default settings are used here, but see https://github.com/linneas/condetri for information and parameters.

# Reset file suffixes to call adapter trimmed files.
fq1=_adaptertrimmed1.fastq
fq2=_adaptertrimmed2.fastq

for samp in $samplist
do
echo "condetri trimming of $samp at ${now}"
perl $condetri -fastq1=${samp}${fq1} -fastq2=${samp}${fq2} -sc=33 -prefix=Qualtrim_${samp}

# Again, optional QC step.
#echo "Performing QC on quality trimmed ${samp} at ${now}"
#fastqc Qualtrim_${samp}_trim1.fastq
#fastqc Qualtrim_${samp}_trim2.fastq
done 

################
# 5) Read deduplication with condetri using default settings.

for samp in $samplist
do
echo "condetri deduplication beginning for ${samp} at $now"
perl $dedup -fastq1=Qualtrim_${samp}_trim1.fastq -fastq2=Qualtrim_${samp}_trim2.fastq -prefix=Dedup_${samp}
done 

# 6) Final QC step

echo "Performing QC on deduplicated ${samp} at $now"
fastqc Dedup_${samp}_uniq1.fastq 
fastqc Dedup_${samp}_uniq2.fastq
done

###################
# At this point you want to check all QC files - pull these down and aggregate.
cat Qualtrim_${samp}*.stats > AllTrim.stats 
cat Qualtrim_${samp}*.hist > AllDedup.hist 
# Then pull down outputs including fastqc files.
#scp -rp [PORT_ID] [PATH_TO_FILES]*.qc [PATH_TO_NEW_LOCATION]


# Clean up intermediate files. 
echo "Cleaning up intermediate files at ${now}"
rm *_adaptertrimmed*.fastq
rm Qualtrim_*_trim*.fastq
#rm *.stats
#rm *.hist
#rm qc files

# Optional: Zip output fastqs if required
#echo "Compressing output fastqs at ${now}"
#for samp in $samplist
#do
# gzip Dedup_${samp}_uniq1.fastq
# gzip Dedup_${samp}_uniq2.fastq
#done

echo "Raw_read_QC_trim.sh completed at ${now}"

