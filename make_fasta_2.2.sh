chmod a+x  make_fasta_2.2_NeSI_reseq.sh
#call FASTA files from BAM files produced by fastq_pipeline_loop[version].sh
#Uses GATK's HaploCaller to call best VARIANTS at any depth
#Nullifies any bases (ref or var) at or less than depth specified by the 'mincoverage' variable

######## variables to change ############
reffile=Kaki_mtgenome_DNA1914_canu_polished_insert.fasta
refdir=/nesi/project/uoo02327/natalie/reference/
ref=$refdir$reffile
#directory that contains the FASTA reference file
grpname=reseq_modern_mt
#name of this group of samples
projectdir=/nesi/nobackup/uoo02327/natalie/reseq_mito_out/
#project directory, that contains data processed by the fastq_pipeline_loop[version].sh script
ploidy=1
#the ploidy of the data (eg, 1 for mitochondria, 2 for autosomal)
dna="modern"
#which BAM file to use
mincoverage="5"
#MT
mtchr="Kaki1mtgenome"
#the MINIMUM coverage allowed to call a base for these samples (eg, if set to 1, all bases covered by 1 or more reads will be in FASTA)
#below: a list of sample names (in any order) as they were processed by the fastq_pipeline_loop[version].sh
#these sample names currently should be in double quotes (around all) and separated by whitespace (spaces or tabs)
samplist="H01383-L1_S41_L007 H01384-L1_S42_L007 H01385-L1_S43_L007 H01386-L1_S1_L001 H01387-L1_S2_L001 H01388-L1_S3_L001 \
H01389-L1_S4_L003 H01390-L1_S5_L003 H01391-L1_S6_L003 H01392-L1_S7_L005 H01393-L1_S8_L005 H01394-L1_S9_L005 H01407-L1_S35_L006 \
H01408-L1_S36_L006 H01409-L1_S37_L006 H01410-L1_S38_L006 H01411-L1_S39_L006 H01412-L1_S40_L006 H01413-L1_S10_L007 H01414-L1_S11_L007 \
H01415-L1_S12_L007 H01416-L1_S1_L001 H01417-L1_S2_L001 H01418-L1_S3_L001"

###########################################################################################################################################################################

#variables that shouldn't need changing
bamlist=$grpname.bam.list
#oneunder=$((mincoverage - 1))
gencalldir="call_genotypes/"
grpcalldir="call_group_genotypes2.0/"

#software versions to use

module load SAMtools/1.9-GCC-7.4.0 snpEff/4.2 GATK/3.8-1 picard/2.1.0 BCFtools/1.9-GCC-7.4.0
SnpSift='/opt/nesi/mahuika/snpEff/4.2/SnpSift.jar'
picard='/opt/nesi/mahuika/picard/2.1.0/picard.jar'
gatk='/opt/nesi/mahuika/GATK/3.8-1/GenomeAnalysisTK.jar'
Rscript='/opt/nesi/mahuika/R/3.6.1-gimkl-2018b/bin/Rscript'
###########################################################################################################################################################################

#enter the project directory

if [[ ! -f /nesi/project/uoo02327/natalie/scripts/filter_cover_for_vcfs.R ]] ; then
    echo "Could not find filter_cover_for_vcfs.R in main processing directory, aborting."
    exit
fi

cd $projectdir

#specify the name of the BAM file to call from
if [ $dna = "ancient" ]; then
echo "You have specified ANCIENT DNA"
bamprefix="_merged_rescale_rdgrps.bam"
elif [ $dna = "modern" ]; then
echo "You have specified MODERN DNA"
bamprefix="_remdup.bam"
fi

#create dictionary file of the reference if it doesn't exist
dictfile=${reffile%.*}.dict
if [ ! -e "$refdir${reffile%.*}.dict" ]; then
echo "Dictionary file of reference does not exist: creating dictionary"
java -jar $picard CreateSequenceDictionary R=$ref O=$refdir${reffile%.*}.dict
else
echo "Dictionary file found"
fi

#index the reference if it is not indexed already
idxfile=${reffile%.*}.fai
if [ ! -e "$refdir$reffile.fai" ]; then
echo "Reference is not indexed: indexing with SAMTOOLS"
samtools faidx $ref
else
echo "Index for reference file found"
fi

#call variants across entire group
if [ ! -e "${projectdir}${grpcalldir}${grpname}_samples.vcf" ]; then
echo "No VCF of called variants found for this group! Please run variantcall_group_vcf_from_bam script first"
else
echo "Using the called variants from ${grpname}_samples.vcf"
fi

###########################################################################################################################################################################


#begin processing each sample
for samp in $samplist
do
#enter directory
echo "starting FASTA creation for sample $samp"

cd ${samp}_process
mkdir $gencalldir
cd $gencalldir

if [ -e "vcf_list.list" ]; then
echo "Removing an existing list of VCF coverage files"
rm vcf_list.list
fi

#get the VCF for the BEST variant calls
echo "creating a VCF of $samp from the group VCF"

java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${projectdir}${grpcalldir}${grpname}_samples.vcf \
-sn $samp \
-env \
-o ${samp}_haploAll.vcf


#make VCF
echo "looking for the BP_RESOLUTION VCF for $samp"

if [ -e "${projectdir}${grpcalldir}${samp}.g.forbam.vcf" ]; then
echo "found BP_RESOLUTION VCF for $samp"
cp ${projectdir}${grpcalldir}${samp}.g.forbam.vcf .
BPvcf=${samp}.g.forbam.vcf
else 
echo "Cannot find BP_RESOLUTION VCF for $samp! Creating a VCF of the entire reference for $samp"
java -jar $gatk \
-T HaplotypeCaller \
-R $ref \
-I ../${samp}${bamprefix} \
-ploidy 1 \
-ERC BP_RESOLUTION \
-o $samp.v.vcf
BPvcf=$samp.v.vcf
fi
#

#make a hack vcf from the ERC vcf
echo "selecting low coverage variants from the VCF for $samp"

#get the VCF
vcfheadl=$(awk '/#CHROM/{ print NR; exit }' $BPvcf)

#get coverage
java -jar $gatk \
 -T DepthOfCoverage \
 -R $ref \
 -o ${samp}_coverage_gatk_table \
 --outputFormat table \
 -I ../../call_group_genotypes2.0/${samp}.gatk.bam

#get the vcf without the header
tail -n +$vcfheadl $BPvcf > $samp.nohead.vcf

#run the VCFs

for cov in $mincoverage
do
#run the script to create a min-cov vcf
echo "running the min.cov script for coverage $cov"
$Rscript /nesi/project/uoo02327/natalie/scripts/filter_cover_for_vcfs.R ${projectdir}${samp}_process/${samp}_coverage_files/${samp}_coverage.txt $samp.nohead.vcf $cov

#
BPvcf=${samp}.g.forbam.vcf

#recreate the vcf
head -n $vcfheadl $BPvcf > $samp.new$cov.vcf
cat $samp.new$cov.vcf gens.$cov.vcf > $samp.new$cov.full.vcf

#remove the NON_REFs from the vcfs
vcffile=$samp.new$cov.full.vcf
echo "altering $vcffile"
sed -i 's/,<NON_REF>//g' $vcffile
sed -i 's/A\t<NON_REF>/N\tA/g' $vcffile
sed -i 's/C\t<NON_REF>/N\tC/g' $vcffile
sed -i 's/G\t<NON_REF>/N\tG/g' $vcffile
sed -i 's/T\t<NON_REF>/N\tT/g' $vcffile


#make fasta for each coverage
echo "creating a FASTA of bases above $cov reads for $samp"

java -Xmx2g -jar $gatk \
-T FastaAlternateReferenceMaker \
-R $ref \
--variant ${samp}_haploAll.vcf \
-snpmask $samp.new$cov.full.vcf \
-snpmaskPriority \
-L $mtchr \
-o ${samp}.min${cov}.filt.fasta
done 

cd ..
cd ..
done


