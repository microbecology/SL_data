#!/bin/bash
#$ -S /bin/bash


source ~/miniconda3/bin/activate
conda activate anvio-6

# directories ## create approprite directories in appropriate path.

WORK=/isilon/ottawa-rdc/users/shared/chenw_lab/galen/sean_li
RAW=/isilon/ottawa-rdc/users/shared/chenw_lab/galen_SeanLi/RawData
origin=$RAW
output=$WORK/output
MAPPING=$output/mapping_anvio




####################################################################
### removal of adapter using sickle
### 
###################################################################


output_folder="sickle"

if [[ ! -d $output/$output_folder ]]; then 
  mkdir -p $output/$output_folder
fi

cd $output/$output_folder

for i in `awk '{print $1}' sample_list`;
do
echo $i
file1="${i}_1.zcat.fastq.gz"
file2="${i}_2.zcat.fastq.gz"
sickle pe -t sanger -f $origin/$file1  -r $origin/$file2 \
                                  -o $output/$output_folder/$i.pair1.fastq \
                                  -p $output/$output_folder/$i.pair2.fastq \
                                  -s $output/$output_folder/$i.single.fastq -q 20 
done


####################################################################
### co-assembly megahit
### 
####################################################################

# input_folder=$WORK/GRDI_QC
# output_folder=$output_dir/megahit_coassembly_anvio

# F=$(echo $(ls $input_folder/*QUALITY_PASSED_R1*)  | sed "s/ /,/g")
# R=$(echo $(ls $input_folder/*QUALITY_PASSED_R2*)  | sed "s/ /,/g")

###run megahit (co-assembly)

# megahit -1 $F -2 $R -o $output_dir/megahit_coassembly_anvio/ --min-contig-len 1000 -m 0.85 --continue -t 16

####################################################################
### mapping with bowtie2
### 
####################################################################

# cd $output_dir/
# mkdir mapping_anvio
# MAPPING=$output_dir/mapping_anvio




# cd $MAPPING

# bowtie2-build $output_dir/megahit_coassembly_anvio/contigs.fa $MAPPING/anvio_contigs --threads 8

# cd $WORK/GRDI_QC
# for sample in `awk '{print $1}' sample_list_anvio.txt`;
# do
# file_ext1="-QUALITY_PASSED_R1.fastq"
# file_ext2="-QUALITY_PASSED_R2.fastq"
# file_ext3=".sam"
# R1="$sample-$file_ext1"
# R2="$sample-$file_ext2"
# sam="$sample-$file_ext3" 
# bowtie2 -x $MAPPING/anvio_contigs -1 $sample$file_ext1 -2 $sample$file_ext2 -S $MAPPING/$sample$file_ext3 --threads 8
# done


### convert sam to bam

# cd $WORK/GRDI_QC

# for sample in `awk '{print $1}' sample_list_anvio.txt`;
# do
# sam=".sam"
# bam=".bam"
# samtools view -S -b $MAPPING/$sample$sam > $MAPPING/$sample$bam
# done

## bam to anvio_bam profiling

# for sample in `awk '{print $1}' sample_list_anvio.txt`;
# do
# anvi-init-bam $MAPPING/${sample}.bam -o $MAPPING/${sample}_anvi.bam
# done


####################################################################
### loading into anvio
### 
####################################################################




#anvi-script-reformat-fasta $output_dir/megahit_coassembly_anvio/final.contigs.fa -o contigs-fixed.fa -l 0 --simplify-names
# mv contigs-fixed.fa contigs.fa

###create anvio database

# anvi-gen-contigs-database -f $output_dir/megahit_coassembly_anvio/contigs.fa -o $output_dir/anvio_anvio/contigs.db -n "GRDI-metagenomics"


# anvi-run-hmms -c $output_dir/anvio_anvio/contigs.db

# anvi-display-contigs-stats $output_dir/anvio_anvio/contigs.db

# anvi-setup-ncbi-cogs $output_dir/anvio_anvio/contigs.db --num-threads 8

# anvi-get-sequences-for-gene-calls -c $output_dir/anvio_anvio/contigs.db \
                                    --get-aa-sequences \
                                    -o $output_dir/anvio_anvio/amino-acid-sequences.fa

sudo anvi-run-ncbi-cogs -c /media/wenRAID/galen/GRDI/output/anvio_sicke/contigs.db -T 20 --cog-data-dir /media/wenRAID/galen/GRDI/output/anvio_anvio/COG

anvi-display-contigs-stats $output_dir/anvio_sicke/contigs.db  --report-as-text  --output-file $output_dir/anvio_sicke/contigs_post-hmm-cogs.txt

####################################################################
###amino-acid seq from eggnog-mapper on SOILMICROBIOME SERVER
### 
####################################################################

#### 

#directories
# ROOT=/media/wenRAID
# GRDI=$ROOT/galen/GRDI
# SOFTWARE=$ROOT/software
# origin_dir=$GRDI/raw_data
# output_dir=$GRDI/output
# input_dir=$GRDI/input
# MAPPING=$output_dir/mapping


# /home/linuxbrew/.linuxbrew/bin/python /media/wenRAID/software/eggnog-mapper/emapper.py -i $output_dir/anvio_anvio/amino-acid-sequences.fa \
#                    --output $output_dir/anvio_anvio/eggNOG/contigs_amino.acid_egg.ann \
#                    --cpu 8 \
#                    -d /media/wenRAID/References/eggnog_DB/hmmdb_levels/NOG_hmm/NOG_hmm.all_hmm


####################################################################
### interproscan on amino-acid seq from anvio
### 
####################################################################


$WORK/interproscan-5.39-77.0/interproscan.sh -i $output_dir/anvio_anvio/amino-acid-sequences.fa -f tsv \
	-d $output_dir/anvio_anvio/interpro/ \
	--tempdir $WORK/temp/ \
	--disable-precalc \
	-appl Pfam,PIRSF,SUPERFAMILY,TIGRFAM \
	--iprlookup \
	--goterms \
	--pathways
