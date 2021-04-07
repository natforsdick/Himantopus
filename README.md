# Himantopus
Scripts related to Forsdick (2020), 'Assessment of the impacts of anthropogenic hybridisation in a threatened non-model bird species through the development of genomic resources with implications for conservation,' PhD thesis research associated with hybridisation between kakī (*Himantopus novaezelandiae*) and poaka (*Himantopus himantopus leucocephalus*).
Thesis permanent link: http://hdl.handle.net/10523/10268

'Raw_read_QC_trim.sh' is a script for performing quality control for short-read whole genome sequencing data. It uses FastQC to assess raw read quality, Trimmomatic for adapter removal, and ConDeTri for quality trimming and read deduplication.

'Streamlined_GBS_pipeline.txt' describes a straightforward analysis pipeline for GBS data including demultiplexing, mapping, variant discovery and filtering, and additional pre-processing prior to analysis with ADMIXTURE and/or adegenet.

'admixture_array_nesi.sh' is a SLURM script to run ADMIXTURE, with ten iterations for each value of K from 1 - 10 processed as an array. 
