chmod a+x fastq_pipeline_loop_modern_WGS.2.5.sh

# Adapted from script by Sophia Cameron-Christie on 20/8/19 by Natalie Forsdick for WGS data
#This script expects a make_coverage_plots Rscript to be in the processing directory
#This script processes multiple paired-end fastq files to merged bam files

#variables
# reference file
reffile=fixedUO_Hhiman_2.0.fasta
# path to reference file
refdir=/nesi/project/uoo02327/natalie/reference/
ref=$refdir$reffile

# path to raw fastq data
datadir=/nesi/nobackup/uoo02327/natalie/data/
home=/nesi/nobackup/uoo02327/natalie/Aus1/map_to_super/
# suffix pattern of the first fastq file
fq1=_fq1.fastq
# suffix pattern of the first fastq file
fq2=_fq2.fastq
# list of samples (separated by spaces, with one pair of quotes around all)
samplist="Trimmed_Aus1_all Trimmed_Aus2_all"

# the name of the target sequence in the FASTA file to look at coverage plots (everything after > and before the first space)

### NEED TO ALTER THIS 
mtchr="contig_1"

#variables that shouldn't need to change
platform="Illumina"
library="Lib1" #currently the library is just set at Lib1 at default, this may be customisable if it is necessary in the future


########################################################################################
#software

## modules are loaded in job script, but as a reference here are the module names to load
module load SAMtools/1.9-GCC-7.4.0
module load picard/2.1.0
module load BWA/0.7.17-gimkl-2017a
module load R/3.5.0-gimkl-2017a

# still need to point to picard file after loading module
picard='/opt/nesi/mahuika/picard/2.1.0/picard.jar'

########################################################################################


#index the reference fasta file
if [ ! -e $refdir$reffile.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
bwa index $ref
else
echo "BWA Index file found"
fi

#Create a Sequence Dictionary if necessary
if [ ! -e $refdir$reffile.dict ]; then
echo "SequenceDictionary of reference does not exist: creating one with picard"
java -jar $picard CreateSequenceDictionary \
R=$ref \
O=$ref.dict
else
echo "SequenceDictionary found"
fi


#####################################################

for samp in $samplist

do

echo "processing data for $samp"

#unzip data if necessary

if [ -e $datadir${samp}${fq1}.gz ] || [ -e $datadir${samp}${fq2}.gz ]; then
echo "Unzipping fastq files"
gunzip $datadir${samp}${fq1}.gz
gunzip $datadir${samp}${fq2}.gz
fi


#get readgroup info
##split the machine/lane info from the fq file
infoline=$(head -n 1 ${datadir}$samp$fq1)
instrument=`echo $infoline | cut -d ':' -f1`
instrument=$(echo $instrument | sed -e 's/ /_/g') # this is to fix issues with spaces in the instrument name
instrumentrun=`echo $infoline | cut -d ':' -f2` 
flowcell=`echo $infoline | cut -d ':' -f3`
lane=`echo $infoline | cut -d ':' -f4`
index=`echo $infoline | cut -d ':' -f10`
#work that info into some 
rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
rgpl="PL:${platform}"
rgpu="PU:${flowcell}.${lane}"
rglb="LB:${samp}_${library}"
rgsm="SM:${samp}"

echo "aligning and sorting data for $samp"

#specify the paired-end reads without the adaptor removals
alignfastq1=${datadir}$samp$fq1
alignfastq2=${datadir}$samp$fq2

#align paired-end files with bwa mem ### HUGH - the original input was from StdIn, not sure what this should be on NeSI
bwa mem -M -t 16 -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rgsm $ref $alignfastq1 $alignfastq2 | java -jar $picard SortSam INPUT=/dev/stdin OUTPUT=${home}${samp}_sort.bam SORT_ORDER=coordinate TMP_DIR=./tmp

#merge the aligned files
alignedbam=${home}${samp}_sort.bam

echo "indexing the bamfile for $samp"

#index the sorted BAM-formatted alignments
samtools index $alignedbam

#get and print the stats from the indexed BAM file (outputs to a text file)
samtools stats $alignedbam >  ${home}${samp}_sort_stats.txt
samtools idxstats $alignedbam > ${home}${samp}_sort_idxstats.txt

## make a plot of the fragment sizes

echo "creating a plot of fragment sizes"

java -jar $picard CollectInsertSizeMetrics \
I=$alignedbam \
O=${home}${samp}_insert_size_metrics.txt \
H=${home}${samp}_insert_size_histogram.pdf \
M=0.5

echo "cleaning the bamfile $samp"

#clean the bam file: Fixes errant MAPQ scores (unmapped reads with scores not zero), soft clips pairs hanging off alignments
java -Xmx50g -jar $picard CleanSam I=$alignedbam O=${home}${samp}_allsortclean.bam CREATE_INDEX=true


#this removes unmapped reads
samtools view -b -F 0x0004 ${home}${samp}_allsortclean.bam > ${home}${samp}_maponly.bam

echo "removing duplicates for $samp"
predup=${home}${samp}_maponly.bam

#remove duplicates## 
java -Xmx50g -jar $picard MarkDuplicates \
I=$predup \
O=${home}${samp}_remdup.bam \
METRICS_FILE=${home}${samp}_remdup_metrics.txt \
CREATE_INDEX=true \
ASSUME_SORTED=true

echo "finished processing $samp"

echo "making plots and stats for $samp"

#create coverage plot

if [ -e /nesi/nobackup/uoo02327/natalie/make_coverage_plots2.0.R ]; then

	echo "Creating coverage plot"

	bamprefix="_remdup.bam"

    mkdir ${home}${samp}_coverage_files

    cd ${home}${samp}_coverage_files

		#create coverage file
		## NEED TO ALTER THIS - pick a chromosome (or several) to map?
		samtools depth -r $mtchr -a ${home}${samp}${bamprefix} > ${home}${samp}_coverage_files/${samp}_coverage.txt

		#run the Rscript to produce a plot

		Rscript /nesi/nobackup/uoo02327/natalie/make_coverage_plots2.0.R $samp $mtchr

	cd ..

else

echo "Cannot find make_coverage_plots.R: skipping coverage plot step"

fi


cd ..

done

##############################
