chmod a+x fastq_pipeline_loop_historic.sh
#13-10-2016 protocol by Catherine Collins #scripting by Sophia Cameron-Christie
# Modified for NeSI on 27/11/19 by Nat Forsdick

#This script processes multiple paired-end fastq files to merged bam files
#OBLIGATORY CHANGES
#change variables
#list sample names
#OPTIONAL CHANGES
#software versions can be specified if necessary
#software versions have all been updated 11/7/2019

#variables
#name of the reference file
reffile=Kaki_mtgenome_DNA1914_canu_polished_insert.fasta
#directory with the reference file
refdir=/nesi/project/uoo02327/natalie/reference/
ref=$refdir$reffile
#directory with the raw fastq data
datadir=/nesi/project/uoo02327/natalie/data/mito/OG4983_fastq/
#suffix pattern of the first fastq file
fq1=_R1.fastq
#suffix pattern of the first fastq file
fq2=_R2.fastq
#list of samples (separated by spaces, with one pair of quotes around all)
samplist="OG4983_MS10992 OG4983_MS11004 OG4983_MS11011"
#excluding samples with too low coverage for successful mapping (MS11008 MS11012 MS11013)
#the name of the target mt sequence in the FASTA file (everything after > and before the first space)
mtchr="Kaki1mtgenome"

#variables that shouldn't need changing
#software/set variables
module purge
module load SAMtools/1.9-GCC-7.4.0 \
picard/2.1.0 \
BWA/0.7.17-GCC-7.4.0 \
Python/2.7.16-gimkl-2018b \
seqtk/1.3-gimkl-2018b \
GSL/2.4-GCCcore-7.4.0 \
AdapterRemoval

mapDamage='/nesi/project/uoo02327/natalie/CondaEnv/mapDamage2/bin/mapDamage'
picard='/opt/nesi/mahuika/picard/2.1.0/picard.jar'
Rscript='/opt/nesi/mahuika/R/3.6.1-gimkl-2018b/bin/Rscript'
export PATH=$PATH:/opt/nesi/mahuika/R/3.6.1-gimkl-2018b/bin/
export R_LIBS=$R_LIBS:/nesi/project/uoo02327/natalie/Rpackages/3.6.1/
covver="2.0"
plotreadver=""
dna="ancient"


#index the reference fasta file
#$bwa index $ref
#index the reference fasta file
if [ ! -e $refdir$reffile.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
bwa index $ref
else
echo "BWA Index file found"
fi


#####################################################

for samp in $samplist
#for samp in MS10560 ###list samples here
do

#mkdir for sample
mkdir /nesi/nobackup/uoo02327/natalie/ancient_mito_out/${samp}_process
cd /nesi/nobackup/uoo02327/natalie/ancient_mito_out/${samp}_process

#unzip data if necessary
if [ -e $datadir${samp}${fq1}.gz ] || [ -e $datadir${samp}${fq2}.gz ]; then
echo "Unzipping fastq files"
gunzip $datadir${samp}${fq1}.gz
gunzip $datadir${samp}${fq2}.gz
fi

#adaptor removal
echo "Running Adapter Removal for "$samp

AdapterRemoval \
--file1 ${datadir}$samp$fq1 \
--file2 ${datadir}$samp$fq2  \
--collapse  \
--trimns  \
--trimqualities  \
--minlength 25  \
--mm 3  \
--minquality 20  \
--basename $samp

#script changed to use .collapsed AdapterRemoval output file
#Find the SA coordinates of the input reads
echo "Finding SA coordinates for "$samp

bwa aln -n 0.03 -o 2 -l 1024 $ref ${samp}.pair1.truncated > ${samp}_p1.sai
bwa aln -n 0.03 -o 2 -l 1024 $ref ${samp}.pair2.truncated > ${samp}_p2.sai
bwa aln -n 0.03 -o 2 -l 1024 $ref ${samp}.collapsed > ${samp}_collapsed.sai

## Repeat the above line for pair2.truncated and collapsed.truncated

#Generate alignments in the SAM format given PAIRED-END reads

bwa sampe $ref ${samp}_p1.sai ${samp}_p2.sai ${samp}.pair1.truncated ${samp}.pair2.truncated > ${samp}.sam

#Generate alignments in the SAM format given SINGLE-END reads

bwa samse $ref ${samp}_collapsed.sai ${samp}.collapsed > ${samp}_collapsed.sam

#get the total reads for calculating endogenous percent

collapsedtotal=$(samtools view -c ${samp}_collapsed.sam)
pairedtotal=$(samtools view -c ${samp}.sam)
echo "$samp $pairedtotal $collapsedtotal" > total_reads.txt

#You will need to do the following for both merged and unmerged reads

#begin loop
for tmpname in ${samp} ${samp}_collapsed
do

echo "Producing filtered bam files for "$samp

#output the alignments as bam, ignoring alignment with quality score lower than 20
samtools view -b -S -q 20 ${tmpname}.sam > ${tmpname}.bam

#sort the alignments by coordinates
samtools sort ${tmpname}.bam -o ${tmpname}_sort.bam

#index the sorted BAM-formatted alignments
samtools index ${tmpname}_sort.bam

#this removes unmapped reads, they should have been removed with q -20 in the first samtools view step
samtools view -b -F 0x0004 ${tmpname}_sort.bam $mtchr > ${tmpname}_maponly.bam 

# Anna added this to remove duplicates, as samtools will remove duplicates from collapsed and uncollapsed reads in a single step these days
samtools rmdup -S --reference $ref ${tmpname}_maponly.bam ${tmpname}_remdup.bam

samtools index ${tmpname}_remdup.bam

#get and print the stats from the indexed BAM file (outputs to a text file)
samtools idxstats ${tmpname}_remdup.bam > ${tmpname}_sort_stats.txt

#end loop
done
#remove duplicates
echo "Removing duplicates for "$samp

java -Xmx8g -jar $picard MarkDuplicates \
I=${samp}_maponly.bam \
remove_duplicates=yes \
O=${samp}_remdup.bam \
METRICS_FILE=${samp}.txt \
ASSUME_SORTED=true

#PALEOMIX remove duplicates from collapsed reads
# Modify for NeSI
echo "Running Paleomix for "$samp
paleomix rmdup_collapsed \
--remove-duplicates \
${samp}_collapsed_maponly.bam > ${samp}_collapsed_remdup.bam

echo "Running mapDamage for "$samp

for tmpname in ${samp} ${samp}_collapsed
do
$mapDamage -i ${tmpname}_remdup.bam -r $ref --rescale 
done

#merge the two bam files after mapDamage
if [ -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "found both rescaled bam files, merging"
samtools merge ${samp}_merged_rescale.bam results_${samp}_remdup/${samp}_remdup.rescaled.bam results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam
bamForReadGroup=${samp}_merged_rescale.bam
#samtools -f merge ${samp}_merged.bam ${samp}_remdup.bam ${samp}_collapsed_remdup.bam
#IF ONLY THE COLLAPSED FILE EXISTS
elif [ ! -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "found the collapsed rescaled bam files, merging"
bamForReadGroup=results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam
#IF ONLY THE PAIRED-END FILE EXISTS
elif [ -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ ! -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "found the paired-end rescaled bam files, merging"
bamForReadGroup=results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam
#IF NEITHER FILE EXISTS
elif [ ! -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ ! -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "cannot find the rescaled bams from MapDamage, check MapDamage has run on either the paired-end or collapsed files"
fi

#add the read groups (required for GATK)

java -Xmx8g -jar $picard AddOrReplaceReadGroups \
I=$bamForReadGroup \
O=${samp}_merged_rescale_rdgrps.bam \
RGID=1 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=$samp

samtools index ${samp}_merged_rescale_rdgrps.bam

echo "finished processing $samp"

echo "making plots and stats for $samp"


echo $samp >> ../contamination_stats.txt
echo "paired-end" >> ../contamination_stats.txt
cat ${samp}_sort_stats.txt >> ../contamination_stats.txt

echo "collapsed" >> ../contamination_stats.txt
cat ${samp}_collapsed_sort_stats.txt >> ../contamination_stats.txt

#create coverage plot

if [ -e /nesi/project/uoo02327/natalie/scripts/make_coverage_plots${covver}.R ]; then

        echo "Creating coverage plot for "$samp

        bamprefix="_merged_rescale_rdgrps.bam"

        mkdir /nesi/nobackup/uoo02327/natalie/ancient_mito_out/${samp}_process/${samp}_coverage_files

        cd /nesi/nobackup/uoo02327/natalie/ancient_mito_out/${samp}_process/${samp}_coverage_files

                #create coverage file
                samtools depth -r $mtchr -a /nesi/nobackup/uoo02327/natalie/ancient_mito_out/${samp}_process/${samp}${bamprefix} > ${samp}_coverage.txt

                #run the Rscript to produce a plot
                $Rscript /nesi/project/uoo02327/natalie/scripts/make_coverage_plots${covver}.R $samp $mtchr # will have to run Rscript separately due to issue with GCC versioning

        cd ..

else

echo "Cannot find make_coverage_plots${covver}.R: skipping coverage plot step"

fi

#create readlength plot

if [ -e /nesi/project/uoo02327/natalie/scripts/plot_readlengths${plotreadver}.R ]; then

        echo "Creating readlength plot for "$samp
                
        cd /nesi/nobackup/uoo02327/natalie/ancient_mito_out/${samp}_process/results_${samp}_collapsed_remdup/

                #run the Rscript to produce a plot
                $Rscript /nesi/project/uoo02327/natalie/scripts/plot_readlengths${plotreadver}.R $samp $PWD # will have to run Rscript separately due to issue with GCC versioning

        cd ..

else

echo "Cannot find plot_readlengths${plotreadver}.R: skipping read length plot step"

fi


cd ..

done


