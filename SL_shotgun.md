SL\_shotgun
================

## loading software (anvio6)

    source ~/miniconda3/bin/activate
    conda activate anvio-6

## create approprite directories in appropriate path.

    GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen/
    WORK=/isilon/ottawa-rdc/users/shared/chenw_lab/galen/sean_li
    RAW=/isilon/ottawa-rdc/users/shared/chenw_lab/galen_SeanLi/RawData
    output=$WORK/output
    TRIM=$output/trimmomatic
    MAPPING=$output/mapping_trimmomatic
    MEGAHIT=$output/megahit_coassembly_trimmomatic
    ANVIO=$output/anvio_trimmomatic

## removal of adapter using trimmomatic

    mkdir $TRIM


    cd 
    for i in `awk '{print $1}' $TRIM/sample_list`;
    do
    echo $i
    file1="${i}_1.zcat.fastq.gz"
    file2="${i}_2.zcat.fastq.gz"
    java -jar $GALEN/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 -trimlog $i.log $RAW/$file1  $RAW/$file2 $TRIM/$i.pair1.fastq.gz $TRIM/$i.unpair1.fastq.gz $TRIM/$i.pair2.fastq.gz $TRIM/$i.unpair2.fastq.gz ILLUMINACLIP:$GALEN/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:6 
    done

    mkdir fastqc_out
    fastqc -t 10 $TRIM/*.pair* $TRIM/fastqc_out

## co-assembly megahit

    F=$(echo $(ls $WORK/output/trimmomatic/*.pair1.fastq.gz)  | sed "s/ /,/g")
    R=$(echo $(ls $WORK/output/trimmomatic/*.pair2.fastq.gz)  | sed "s/ /,/g")


    megahit -1 $F -2 $R -o $MEGAHIT -t 18 --min-contig-len 1000 -m 0.85

## mapping with bowtie2

    cd $output/
    mkdir mapping_trimmomatic
    MAPPING=$output/mapping_trimmomatic


    cd $MAPPING

    bowtie2-build $MEGAHIT/contigs.fa $MAPPING/trimmomatic_contigs --threads 20

    cd 
    for i in `awk '{print $1}' $TRIM/sample_list`;
    do
    echo $i
    file1="${i}_pair1.fastq.gz"
    file2="${i}_pair2.fastq.gz"
    sam="${i}.sam"
    bowtie2 -x $MAPPING/anvio_contigs -1 $file1 -2 $file2 -S $MAPPING/$sample$file_ext3 --threads 20
    done

convert sam to bam

    for sample in `awk '{print $1}' sample_list`;
    do
    sam=".sam"
    bam=".bam"
    samtools view -S -b $MAPPING/$sample$sam > $MAPPING/$sample$bam
    done

bam to anvio\_bam profiling

    for sample in `awk '{print $1}' sample_list`;
    do
    anvi-init-bam $MAPPING/${sample}.bam -o $MAPPING/${sample}_anvi.bam
    done

## loading into anvio

    anvi-script-reformat-fasta $MEGAHIT/final.contigs.fa -o $MEGAHIT/contigs-fixed.fa -l 1000 --simplify-names
    mv $MEGAHIT/contigs-fixed.fa $MEGAHIT/contigs.fa

create anvio database

    mkdir $output/anvio_trimmomatic
    anvi-gen-contigs-database -f $MEGAHIT/contigs.fa -o $ANVIO/contigs.db -n "SeanLi-metagenomics"

    anvi-run-hmms -c $ANVIO/contigs.db

    anvi-display-contigs-stats $ANVIO/contigs.db --report-as-text  --output-file $ANVIO/contigs_post-hmm-cogs.txt

    anvi-setup-ncbi-cogs $ANVIO/contigs.db --num-threads 8

    anvi-run-ncbi-cogs -c $ANVIO/contigs.db --cog-data-dir $GALEN/COG --num-threads 20

    anvi-get-sequences-for-gene-calls -c $ANVIO/contigs.db \
                                        --get-aa-sequences \
                                        -o $ANVIO/amino-acid-sequences.fa

## interproscan on amino-acid seq from anvio


    $WORK/interproscan-5.39-77.0/interproscan.sh -i $ANVIO/amino-acid-sequences.fa -f tsv \
        -d $ANVIO/interpro/ \
        --tempdir $WORK/temp/ \
        --disable-precalc \
        -appl Pfam,PIRSF,SUPERFAMILY,TIGRFAM \
        --iprlookup \
        --goterms \
        --pathways

create new folder to house concatenate of all smaller fxnal annotation output (the script dont like a folder with too much clutter, im guessing)

    mkdir $ANVIO/interpro/iprs_output

    cat $ANVIO/interpro/*.tsv > $ANVIO/interpro/iprs_output/all_iprs.tsv

script to clean up and allow import to anvio create iprs2anvio.sh file found here: <https://github.com/xvazquezc/stuff/blob/master/iprs2anvio.sh>

    cd $ANVIO/interpro/iprs_output
    /isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.39-77.0/iprs2anvio.sh -i $ANVIO/interpro/iprs_output/all_iprs.tsv -o $ANVIO/interpro/iprs_output/all_iprs -g -p -r


    cat $ANVIO/interpro/*.tsv > $ANVIO/interpro/iprs_output/all_iprs.tsv

importing functional annotation to anvio

    anvi-import-functions -c $ANVIO/contigs.db -i $ANVIO/interpro/iprs_output/all_iprs_iprs2anvio.tsv

## TAXONOMY CALLS

centrifuge on amino-acid seq from anvio

    mkdir $ANVIO/centrifuge

    cd $ANVIO/centrifuge

    centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v $ANVIO/amino-acid-sequences.fa -S $ANVIO/centrifuge/centrifuge_hits.tsv

Make sure there is two files in the work directory ($ANVIO/centrifuge/)

    anvi-import-taxonomy-for-genes -c $ANVIO/contigs.db -i $ANVIO/centrifuge/centrifuge_report.tsv $ANVIO/centrifuge/centrifuge_hits.tsv -p centrifuge

## Profiling cont'd very long. bam --&gt; anvio

PROFILE DONT LIKE SAMPLE NAME WITH "-"

create bash file running


    file_ext="_anvi.bam"
    for sample in `awk '{print $1}' $ANVIO/profile/bam_aa`;
    do
    echo "$sample$file_ext"
    anvi-profile -i $MAPPING/"$sample$file_ext" -c $ANVIO/contigs.db --output-dir $ANVIO/profile/$sample --sample-name $sample -T 10
    done

merge profile (very long...)

    anvi-merge -p $ANVIO/profile/*/PROFILE.db -o $ANVIO/profile_merged -c $ANVIO/contigs.db -S sean_li

## BINNING! Finally

Concoct

    anvi-cluster-contigs -p $ANVIO/profile_merged/PROFILE.db -c $ANVIO/contigs.db -C metagenomic_anvio_concoct --driver concoct -T 50 --just-do-it

metabat2

    anvi-cluster-contigs -p $ANVIO/profile_merged/PROFILE.db -c $ANVIO/contigs.db -C metagenomic_anvio_metabat2 --driver metabat2 -T 50 --just-do-it

maxbin2

    anvi-cluster-contigs -p $ANVIO/profile_merged/PROFILE.db -c $ANVIO/contigs.db -C metagenomic_anvio_maxbin2 --driver metabat2 -T 50 --just-do-it

binsanity

    anvi-cluster-contigs -p $ANVIO/profile_merged/PROFILE.db -c $ANVIO/contigs.db -C metagenomic_binsanity --driver metabat2 -T 50 --just-do-it

dastool

    anvi-cluster-contigs -p $ANVIO/profile_merged/PROFILE.db -c $ANVIO/contigs.db -S metagenomic_anvio_concoct,metagenomic_anvio_metabat2,metagenomic_anvio_maxbin2,metagenomic_binsanity --search-engine diamond -C metagenomic_anvio_dastool -T 50 --just-do-it

checkm

    checkm lineage_wf $ANVIO/sample_summary_1_dastool/bin_by_bin/checkm/ $ANVIO/sample_summary_1_dastool/checkm -t 20 --tmp $ANVIO/sample_summary_1_dastool/tmp -x fa --tmpdir $WORK/temp

## GTDB-tk (taxonomical classification)

    GTDBTK_DATA_PATH="/isilon/ottawa-rdc/users/shared/chenw_lab/galen/database/gtdb/"
    gtdbtk de_novo_wf --genome_dir $ANVIO/sample_summary_1_dastool/bin_by_bin/all_bins --bacteria -x fa --outgroup_taxon p__Patescibacteria --out_dir $ANVIO/sample_summary_1_dastool/gtdb --cpus 20


    gtdbtk ani_rep --genome_dir $ANVIO/sample_summary_1_dastool/bin_by_bin/all_bins --out_dir $ANVIO/sample_summary_1_dastool/gtdb -x fa --cpus 20

    gtdbtk classify_wf --genome_dir $ANVIO/sample_summary_1_dastool/bin_by_bin/all_bins --out_dir $ANVIO/sample_summary_1_dastool/gtdb --outgroup_taxon p__Patescibacteria -x fa --cpus 20

## DRAM (functional profiling)

    cd $ANVIO/sample_summary_1_dastool/bin_by_bin/all_bins/

    DRAM.py annotate -i '*contigs.fa' -o $ANNOTATION --threads 50

Distilling (summarizing all data into a very nice graph)

    DRAM.py distill -i $ANNOTATION/annotations.tsv -o $ANNOTATION/genome_summaries --trna_path $ANNOTATION/trnas.tsv --rrna_path $ANNOTATION/rrnas.tsv
