#!/usr/bin/env bash

# Created by Denise Martini, 2018.
# Modified by Natalie Forsdick, 2018.

# Script to loop over sequence files and map against a reference genome.
# Specifically used for stilt GBS data.

# This script takes two arguments: 1st: a text file with list of fq files to map, 
# and 2nd: the reference genome.
# NB: the arguments have to be in the right order

# like this: 
# mapping_loop.sh fq_list /path/to/ref/ref_genome.fasta
# First, check if the reference genome has been indexed:


if [ ! -f $2\.amb ]
	then
	bwa-0.7.12-r1044 index -a bwtsw $2
fi

# bwa mem -t "${bwathr}" ../refgenome/"${refgen}" {}.fastq ">" {}.sam

# Mapping steps:

cat $1 | while read fq
do
	echo $fq
	base=$(echo $fq | cut -f 1 -d '.')
	bwa-0.7.12-r1044 mem -t 4 $2 $fq > tmp_$base\.sam
	samtools-1.2 view -b -S -h tmp_$base\.sam > tmp_$base\.bam
	samtools-1.2 sort tmp_$base\.bam -o sorted_$base\.bam
	samtools-1.2 index sorted_$base\.bam


done

# Clean up steps:

rm tmp_*

mv *.bam ../mapped
mv *.bai ../mapped
