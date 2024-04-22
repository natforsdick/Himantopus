chmod a+x variantcall_group_vcf_from_bam2.5_NeSI.sh
#call VCF files from BAM files produced by fastq_pipeline_loop[version].sh
#Uses GATK's HaploCaller to call all variants, and combine into a group VCF 

# Protocol by Sophia Cameron-Christie, altered by Nat Forsdick for use on NeSI, 3/12/19.
#variables to change (description below each)

reffile=Kaki_mtgenome_DNA1914_canu_polished_insert.fasta
#directory with the reference file
refdir=/nesi/project/uoo02327/natalie/reference/
ref=$refdir$reffile
#directory that contains the FASTA reference file
grpname=stilt_modern_mt
#name of this group of samples
projectdir=/nesi/nobackup/uoo02327/natalie/ancienct_mito_out/
#project directory, that contains data processed by the fastq_pipeline_loop[version].sh script
ploidy=1
#the ploidy of the data (eg, 1 for mitochondria, 2 for autosomal)
callingtype="split"
#whether to call the variants as "split" or "grouped". "split" is the method recommended by GATK, but the "grouped" method takes significant amounts of time but make work well for aDNA.
#the "grouped" method may not be recommended if:
#if sample number is about 20 this method is not recommended
#if sequence target is anything larger than the mitochondria this method may not be recommended
#if the depth of sequencing is high (>100 average, as in modern mtDNA) this method may not be recommended
filtertype="hard"
#how to filter the variants; either "VQSR" or "hard" 
#for genome-wide data, the recommendation is Variant Quality Score Recalibration (VQSR) but this is not recommended for when only mitochondrial data is available
#because the number of variants are too few; therefore "hard" filtering is used for these variants
#hard filters SHOULD BE EXAMINED AND ADAPTED according to the dataset; further discussion can be found on GATK's site: software.broadinstitute.org/gatk/guide/article?id=6925
dna="ancient"
#the prefix of the .bam file that the coverage should be called from: usuall "_remdup.bam" for modern DNA and "_merged_rescale_rdgrps.bam" for ancient
MQthreshold="30.0"
#the quality threshold for filtering - should be specific to the group # set to 40 initially, can reduce this as required.
call="yes"
#either yes or [anything else]. If "yes", variants will be called from sample BAMs. If it is not yes, existing sample VCFs will be used.

#below: a list of sample names (in any order) as they were processed by the fastq_pipeline_loop[version].sh
#these sample names currently should be in double quotes (around all) and separated by whitespace (spaces or tabs)
samplist="MS10986 MS10987 MS10988 MS10990 MS10991 MS10992 MS10994" 

# First you need to rename all your output files from the previous step to match the MS numbers
# eg (test with 'echo', then replace with 'mv'): 
#for f in fgh*; do echo "$f" $(echo "$f" | sed 's/^fgh/jkl/g'); done

###########################################################################################################################################################################

#variables that shouldn't need changing
ref=$refdir$reffile
grpcalldir="call_group_genotypes2.0"

#software versions to use
module purge
module load SAMtools/1.9-GCC-7.4.0 picard/2.1.0 BCFtools/1.9-GCC-7.4.0 snpEff/4.2 GATK/3.8-1
#samtools=/usr/local/samtools-1.7/samtools
gatk='/opt/nesi/mahuika/GATK/3.8-1/GenomeAnalysisTK.jar'
#module load picard/2.18.0

vcfutils='/opt/nesi/mahuika/BCFtools/1.9-GCC-7.4.0/bin/vcfutils.pl'
#samtools=/usr/local/samtools-1.3.1/bin/samtools
SnpSift='/opt/nesi/mahuika/snpEff/4.2/SnpSift.jar'
picard='/opt/nesi/mahuika/picard/2.1.0/picard.jar'
#gatk=/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar


###########################################################################################################################################################################

# "13" error in reference file solved with dos2unix
#dos2unix $ref
#sed -i 's/-/N/g' $ref
#sed -i 's/n/N/g' $ref

#specify the name of the BAM file to call from
if [ $dna = "ancient" ]; then
echo "You have specified ANCIENT DNA"
bamprefix="_merged_rescale_rdgrps.bam"
elif [ $dna = "modern" ]; then
echo "You have specified MODERN DNA"
bamprefix="_remdup.bam"
fi

if [ ! -d "$grpcalldir" ]; then
echo "creating new directory for group variant calling"
mkdir $grpcalldir
fi

idxfile=${reffile%.*}.fai
if [ ! -e "$refdir$reffile.fai" ]; then
echo "Reference is not indexed: indexing with SAMTOOLS"
samtools faidx $ref
else
echo "Index for reference file found"
fi

dictfile=${reffile%.*}.dict
if [ ! -e "$refdir${reffile%.*}.dict" ]; then
echo "Dictionary file of reference does not exist: creating dictionary"
java -Xmx4g -jar $picard CreateSequenceDictionary R=$ref O=$refdir${reffile%.*}.dict
else
echo "Dictionary file found"
fi

cd $grpcalldir

if [ $callingtype = "split" ]; then

#begin the calling on each sample
if [ $call = "yes" ]; then
echo "variants will be called from individual sample BAMs"
for samp in $samplist
do

# #call sample
#java -jar gatk \
#-T HaplotypeCaller \
#-R $ref \
#-I $projectdir${samp}_process/${samp}${bamprefix} \
#-ploidy 1 \
#--emitRefConfidence BP_RESOLUTION \
#-bamout ${samp}.gatk.bam \
#-o $samp.g.forbam.vcf

java -jar $gatk \
-T HaplotypeCaller \
-R $ref \
-I $projectdir${samp}_process/${samp}${bamprefix} \
-ploidy 1 \
--emitRefConfidence BP_RESOLUTION \
-bamout ${samp}.gatk.bam \
-o $samp.g.forbam.vcf
#
# #

done
else

echo "variants will NOT be called from individual sample BAMs"

fi
#remove the existing list of vcfs

if [ -e "vcf.list" ]; then
echo "removing the existing list of VCFs"
rm vcf.list
fi

#make a list of the gvcfs to compare group genotypes
for samp in $samplist
do
echo $samp.g.forbam.vcf >> vcf.list
done

#joint genotyping

#java -jar
java -jar $gatk \
-T GenotypeGVCFs \
-R $ref \
--variant vcf.list \
-o ${grpname}_unfiltered_variants.vcf

#variant recalibration (not recommended when only using mtDNA vars) /filtering

#filter the SNPs

echo "pulling out SNPs and INDELs"

#select snps
#java -jar $
java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${grpname}_unfiltered_variants.vcf \
-selectType SNP \
-o raw_snps.vcf 

#select the indels
#java -jar $
java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${grpname}_unfiltered_variants.vcf \
-selectType INDEL \
-o raw_indels.vcf 

echo "filtering the SNPs and INDELs"

if [ $filtertype = "hard" ]; then

#filter snps
#java -jar $
java -jar $gatk \
-T VariantFiltration \
-R $ref \
--variant raw_snps.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < $MQthreshold" \
--filterName "gatk_snp_filter" \
-o ${grpname}_filtered_snps.vcf 


#filter indels
#java -jar $
java -jar $gatk \
-T VariantFiltration \
-R $ref \
-V raw_indels.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filterName "gatk_indel_filter" \
-o ${grpname}_filtered_indels.vcf 

#recombine indels + snps

#recombine indels + snps

#java -jar $
java -jar $gatk \
-T CombineVariants \
-R $ref \
--variant:snp ${grpname}_filtered_snps.vcf \
--variant:ind ${grpname}_filtered_indels.vcf \
-genotypeMergeOptions PRIORITIZE \
-priority snp,ind \
-o ${grpname}_samples_include-filtered.vcf

#remove variants that have failed filtering

#java -jar $
java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${grpname}_samples_include-filtered.vcf \
-ef \
-o ${grpname}_samples.vcf


else

echo "currently not set up to do VQSR filtering, please perform manually"

fi

########################################################################################

else

#create bamlist
if [ ! -e "$projectdir$grpname.bam.list" ]; then
echo "Could not find bamlist! Creating one from the samplist"
ls ${projectdir}*_process/*${bamprefix} > $projectdir$grpname.TEST.bam.list
else
echo "Found bamlist for group $grpname"
fi


#call the variants 
if [ ! -e "$projectdir${grpname}_samples.vcf" ]; then
echo "No VCF of called variants for this group! Calling a new VCF with HaplotypeCaller. This make take some time"
#java -jar $
java -jar $gatk \
-T HaplotypeCaller \
-R $ref \
-I ${projectdir}${bamlist} \
-ploidy $ploidy \
-o $projectdir${grpcalldir}${grpname}_samples.vcf
else
echo "VCF of called variants exists for this group. Please delete the existing VCF if you wish to overwrite"
fi
fi

cd ..

