#!/bin/sh -e

#  Remote_blast_multiple.sh
#  
#  Created by NatalieForsdick on 20/03/18.
#
# To check sequence data against the reference database (for contamination etc.).
# First we load Blast as a module.
module load BLAST+

# Specify the samples for checking.
samplist='LIST_OF_SAMPLE_FASTAS'
fasta=.fa

# Compare against the genome database.
for samp in $samplist
do
blastn -db refseq_representative_genomes -query ${samp}${fasta} -out ${samp}.out -remote

done

# If you want to split files so they represent subsamples, that will then be quickly compared against the database:
# Here I was interested in checking the reads that failed to map to the kakī draft reference genome.

samplist='LIST_OF_SAMPLE_FASTAS'

for samp in $samplist
do
echo "Splitting ${samp} fasta"
pyfasta split -n 6 ${samp}_all_PEunmapped2.fa 
pyfasta split -n 6 ${samp}_all_PEunmapped1.fa
done

