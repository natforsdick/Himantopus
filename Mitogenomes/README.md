# Mitogenomes

This directory contains scripts for processing mitogenome data of *Himantopus spp.*, associated with the draft manuscript 'Population decline and conservation management of the critically endangered kakī/black stilt (*Himantopus novaezelandiae*) informed by mitochondrial diversity', Forsdick et al., in progress. These script were originally developed by Sophia Cameron-Christie for use by the Matisoo-Smith lab group in the Department of Anatomy, University of Otago.

For historic samples, I processed all sequencing data with [fastq_pipeline_loop_historic.sh](Mitogenomes/fastq_pipeline_loop_historic.sh). 

For modern samples, I processed all sequencing data with [fastq_pipeline_loop_modern_WGS.2.5.sh](Mitogenomes/fastq_pipeline_loop_modern_WGS.2.5.sh).

These pipelines use the additional R scripts supplied to plot coverage depth across the mitogenomes, and then filter based on this information. 

Variant calling was done using [variantcall_group_vcf_from_bam2.5_NeSI.sh](Mitogenomes/variantcall_group_vcf_from_bam2.5_NeSI.sh).

Finally, output FASTA files of the consensus sequences were produced using [make_fasta_2.2.sh](Mitogenomes/make_fasta_2.2.sh).

