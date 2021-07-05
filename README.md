## Himantopus
Scripts related to Forsdick (2020), 'Assessment of the impacts of anthropogenic hybridisation in a threatened non-model bird species through the development of genomic resources with implications for conservation,' PhD thesis associated with interspecific hybridisation between kakī (*Himantopus novaezelandiae*) and poaka (*Himantopus himantopus leucocephalus*). Thesis permanent link: http://hdl.handle.net/10523/10268.

Associated peer-reviewed publications:
- Galla, Forsdick et al. (2019). Reference Genomes from Distantly Related Species Can Be Used for Discovery of Single Nucleotide Polymorphisms to Inform Conservation Management. _Genes_, 10(1), 9.  https://doi.org/10.3390/genes10010009.
- Forsdick et al. (2021).Genomic sequencing confirms absence of introgression despite past hybridisation between a critically endangered bird and its common congener. _Global Ecology and Conservation_, 28, e01681. https://doi.org/10.1016/j.gecco.2021.e01681. 

# Files
_Raw_read_QC_trim.sh_ is a script for performing quality control for short-read whole genome sequencing data. It uses FastQC to assess raw read quality, Trimmomatic for adapter removal, and ConDeTri for quality trimming and read deduplication.

_Streamlined_GBS_pipeline.md_ describes a straightforward analysis pipeline for GBS data including demultiplexing, mapping, variant discovery and filtering, and additional pre-processing prior to analysis with ADMIXTURE and/or adegenet.

_admixture_array_nesi.sh_ is a SLURM script to run ADMIXTURE, with ten iterations for each value of K from 1 - 10 processed as an array. 

_pophelper-analysis.Rmd_ is an Rmarkdown script to merge and visualise the results of multiple iterations of ADMIXTURE analysis. 

_Forsdick-et-al_Stilts_SNPRelate.Rmd_ includes the methods used to perform population clustering, and estimation of relatedness and Fst using the R package SNPRelate. The results of this analysis are appended as Supplementary File 1 in Forsdick et al. (_In review_).
