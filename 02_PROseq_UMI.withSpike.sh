#!/usr/bin/bash

#bash 01_proseq_preprocessing_UMI.sh /share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/PRO-seq-UMI 20210702 hg19
if [ $# -ne 3 ]; then
    echo -e $0: usage: bash $0 "$1: Path of directory; $2: datainfo and $3:reference genome"
    exit 1
fi


###file_path
all_path=$1
datainfo=$2
reference_info=$3
input_path=${all_path}/Raw_data/${datainfo}
sampleinfo=${input_path}/sampleinfo_${datainfo}.txt
file_path=${all_path}/Alignment/${datainfo}   # Run_Path

raw_dir=${file_path}/00_rawdata
logs_dir=${file_path}/logs
fastqc_dir=${file_path}/01_fastqc
trimmedFastq_dir=${file_path}/02.1_trimmeddata
fastqc_dir=${file_path}/01_fastqc
trimmedFastq_log_dir=${file_path}/logs/trimmeddata
trimmmed_fastp_dir=${file_path}/02.2_fastp_labelUMI
Fastp_log_dir=${file_path}/logs/Fastp_trimmeddata
rmrRNAdata_dir=${file_path}/02.3_rmrRNAdata
rmrRNAdata_log_dir=${file_path}/logs/rRNA
rmrRNAdata_fastqc_dir=${file_path}/02.4_rmrRNAdata_fastqc

align_exp_dir=${file_path}/03_bam
alignexp_log_dir=${file_path}/logs/align_exp
align_spike_dir=${file_path}/03_spikebam
alignspike_log_dir=${file_path}/logs/align_spike
exp_bam_rmdup=${file_path}/03_bam_rmdup
rmdup_exp_log=${file_path}/logs/rmdup_state_exp
spike_bam_rmdup=${file_path}/03_spikebam_rmdup
rmdup_spike_log=${file_path}/logs/rmdup_state_spike


bw_fulllength_dir=${file_path}/04_bw_fulllength_rmdup
bw_singlebase_dir=${file_path}/04_bw_singlebase_rmdup

sampleinfo=${input_path}/sampleinfo_${datainfo}.txt


#reference genome
GENOME_human="/share/home/Blueberry/reference/index/bowtie2/hg19/hg19"
GENOME_mouse="/share/home/Blueberry/reference/index/bowtie2/bowtie2_download/mm10/mm10"
RDNA_human="/share/home/Blueberry/reference/index/bowtie2/rDNA/rDNA"
RDNA_mouse="/share/home/Blueberry/reference/index/bowtie2/mm10_rDNA/mm10_rDNA"

if [[ ${reference_info} == "hg19" ]]
then
    GENOME_EXP=${GENOME_human}
    GENOME_SPIKE=${GENOME_mouse}
    GENOME_RDNA=${RDNA_human}
    exp_info="hg19"
    spike_info="mm10"
else
    GENOME_EXP=${GENOME_mouse}
    GENOME_SPIKE=${GENOME_human}
    GENOME_RDNA=${RDNA_mouse}
    exp_info="mm10"
    spike_info="hg19"
fi

echo -e "\n***************************\nPRO-seq processing at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
echo -e "experimental genome is: ${GENOME_EXP} \nspike-in genome is: ${GENOME_SPIKE} \nrDNA is: ${GENOME_RDNA}"



## Adaptor sequences to clip. Default = Tru-Seq small RNA
ADAPTOR_1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" 
ADAPTOR_2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"


MAPQ=10


#step 1.1
####change file name#####
echo -e "\n***************************\nRenaming files at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${raw_dir}

cat $sampleinfo| while read id;
do
    arr=($id)
    sample1=${arr[0]}
    sample2=${arr[1]}
    echo "input_path:${input_path}"
    if [ ! -s ${raw_dir}/${sample2}_R1.fastq.gz ]
    then
        fq1=$(ls ${input_path}/*/*/*R1*f*q.gz|grep "$sample1")
        fq2=$(ls ${input_path}/*/*/*R2*f*q.gz|grep "$sample1")
        ln -s $fq1 ${raw_dir}/${sample2}_R1.fastq.gz
        ln -s $fq2 ${raw_dir}/${sample2}_R2.fastq.gz
    fi
done


#step 1.2
####fastqc of raw data ####
echo -e "\n***************************\nfastqc of raw data at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${fastqc_dir}

in_path=${raw_dir}
out_path=${fastqc_dir}

nohup_number=0
for FILE in `ls ${in_path}/*.gz`
do 
	sample=$(basename ${FILE/.fastq.gz/})
    if [ ! -s ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})" ]
	then
        echo "Generating file: ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})"..."
	    fastqc $FILE -t 4 -o ${out_path}/ > ${out_path}/${sample}_fastqc.log 2>&1 &
	fi
    nohup_number=`echo $nohup_number+4 | bc`
    if [[ $nohup_number -eq 28 ]]
    then
        echo "waiting..."
        wait
        nohup_number=0
    fi
done



wait
#step 1.3
#check fq
#cd ${fastqc_dir}
#ls *fastqc.zip|sed 's/_R[1-2]_fastqc.zip//g'|sort|uniq -c

### merge reports of fastqc
multiqc ${fastqc_dir}/ -n rawdata_multiqc_${datainfo} -o ${fastqc_dir}/


#step 2.1
### Trimming adapters  (trim_galore)
####remember to change --length for proseq !!!!!! #####
echo -e "\n***************************\nTrimming adapters at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${trimmedFastq_dir}
mkdir -p ${trimmedFastq_log_dir}

nohup_number=0
for fq1 in `ls ${raw_dir}/*R1.fastq.gz`
do
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}
    if [ ! -s ${trimmedFastq_dir}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})" ]
    then
        echo "Generating file: ${trimmedFastq_dir}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})";"
     trim_galore -q 25 --phred33 --length 15 -e 0.1 --stringency 4 --paired -o ${trimmedFastq_dir} $fq1 $fq2 \
    > ${trimmedFastq_log_dir}/"$(basename ${fq1/_R1.fq.gz/_trimmed.log})" 2>&1 &
    fi
    nohup_number=`echo $nohup_number+1 | bc`
    if [[ $nohup_number -eq 28 ]]
    then
        echo "waiting..."
        wait
        nohup_number=0
    fi
done

wait


#step 2.2
### label umi (fastp)
echo -e "\n***************************\nlabeling umi at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"


mkdir -p ${trimmmed_fastp_dir}
mkdir -p ${Fastp_log_dir}

for fq1 in `ls ${trimmedFastq_dir}/*R1_val_1.fq.gz`
do
    fq2=${fq1/R1_val_1.fq.gz/R2_val_2.fq.gz}
    sample=$(basename ${fq1/_R1_val_1.fq.gz/})
    if [ ! -s ${trimmmed_fastp_dir}/${sample}_fastp_R1.fastq.gz ]
    then
    fastp \
                        -i $fq1 \
                        -o ${trimmmed_fastp_dir}/${sample}_fastp_R1.fastq.gz \
                        -I $fq2 \
                        -O ${trimmmed_fastp_dir}/${sample}_fastp_R2.fastq.gz \
                        --adapter_sequence $ADAPTOR_1 \
                        --adapter_sequence_r2 $ADAPTOR_2 \
                        --umi \
                        --umi_loc=per_read \
                        --umi_len=6 \
                        --umi_prefix="UMI" \
                        --html ${Fastp_log_dir}/${sample}_fastp.html \
                        -w 25 \
                        -c \
                        --overlap_len_require 15 2> ${Fastp_log_dir}/${sample}_fastp.log
                    fi

done


wait


#step 2.3
####filtering rRNA reads#####
echo -e "\n***************************\nfiltering rRNA reads at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${rmrRNAdata_dir}
mkdir -p ${rmrRNAdata_log_dir}

#mouse
for fq1 in `ls ${trimmmed_fastp_dir}/*_fastp_R1.fastq.gz`
do
    fq2=${fq1/fastp_R1.fastq.gz/fastp_R2.fastq.gz}
    sample=$(basename ${fq1/_fastp_R1.fastq.gz/})
    if [ ! -s ${rmrRNAdata_log_dir}/${sample}_rRNA_bowtie.log ]
    then
        echo "Generating file: ${rmrRNAdata_log_dir}/${sample}_rRNA_bowtie.log; "
    bowtie2 \
                --fast-local \
                -x ${GENOME_RDNA} \
                -1 $fq1 \
                -2 $fq2 \
                --un-conc-gz  ${rmrRNAdata_dir}/${sample}_rmrRNA.fq.gz \
                --threads 25 2> ${rmrRNAdata_log_dir}/${sample}_rRNA_bowtie.log > /dev/null
    fi

done




wait

#step 2.4
#####rename rmrRNAdata####
for FILE in ${rmrRNAdata_dir}/*.fq.1.gz
do
    if [ ! -s ${FILE/_rmrRNA.fq.1.gz/_R1.fq.gz} ]
    then
        mv "$FILE" ${FILE/_rmrRNA.fq.1.gz/_R1.fq.gz}
    fi
done

for FILE in ${rmrRNAdata_dir}/*.fq.2.gz
do
    if [ ! -s ${FILE/_rmrRNA.fq.2.gz/_R2.fq.gz} ]
    then
        mv "$FILE" ${FILE/_rmrRNA.fq.2.gz/_R2.fq.gz}
    fi
done

wait



#step 2.5
#####qc for trimmed data####
echo -e "\n***************************\nqc for trimmed data at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${rmrRNAdata_fastqc_dir}

nohup_number=0
for FILE in `ls ${rmrRNAdata_dir}/*fq.gz`
do
    sample=$(basename ${FILE/.fq.gz/})
    if [ ! -s ${rmrRNAdata_fastqc_dir}/"$(basename ${FILE/.fq.gz/_fastqc.html})" ]
    then
        echo "Generating file: ${rmrRNAdata_fastqc_dir}/"$(basename ${FILE/.fq.gz/_fastqc.html})"; "
    fastqc $FILE -t 25 -o ${rmrRNAdata_fastqc_dir} >  ${rmrRNAdata_fastqc_dir}/${sample}_fastqc.log 2>&1 &
    fi
    nohup_number=`echo $nohup_number+4 | bc`
    if [[ $nohup_number -eq 28 ]]
    then
        echo "waiting..."
        wait
        nohup_number=0
    fi

done

wait

### merge reports of fastqc
multiqc ${rmrRNAdata_fastqc_dir}/ -n rmrRNA_multiqc_${datainfo} -o ${rmrRNAdata_fastqc_dir}/


#step 3.1
###Aligning to experimental genome#####
echo -e "\n***************************\nAligning to experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p  ${align_exp_dir}
mkdir -p  ${alignexp_log_dir}

for PAIR in $(ls ${rmrRNAdata_dir} | sed 's/_R[1-2].*//' | uniq )
do
if [ ! -s "${align_exp_dir}/${PAIR}_${exp_info}.bam" ]
then
    echo "aligning ${PAIR} to experimental genome"
    (bowtie2 \
    --local \
    --sensitive-local \
    --threads 25 \
    -x "$GENOME_EXP" \
    -1 "${rmrRNAdata_dir}/${PAIR}_R1.fq.gz" \
    -2 "${rmrRNAdata_dir}/${PAIR}_R2.fq.gz" \
    2> ${alignexp_log_dir}/${PAIR}_align.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
    samtools sort -@ 25 -o ${align_exp_dir}/${PAIR}_${exp_info}.bam
    samtools index -@ 25 ${align_exp_dir}/${PAIR}_${exp_info}.bam
fi
done


wait


#step 3.2
### Aligning to spike-in genome to get normalization factors ###
echo -e "\n***************************\nAligning to spike-in genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${align_spike_dir}
mkdir -p ${alignspike_log_dir}

for PAIR in $(ls ${rmrRNAdata_dir} | sed 's/_R[1-2].*//' | uniq )
do
if [ ! -s "${align_spike_dir}/${PAIR}_${spike_info}.bam" ]
    then
    echo "aligning ${PAIR} to spike-in genome"
    (bowtie2 \
    --local \
    --very-sensitive-local \
    --threads 25 \
    --no-unal \
    --no-mixed \
    --no-discordant \
    -x "$GENOME_SPIKE" \
    -1 "${rmrRNAdata_dir}/${PAIR}_R1.fq.gz" \
    -2 "${rmrRNAdata_dir}/${PAIR}_R2.fq.gz" \
    2> ${alignspike_log_dir}/${PAIR}_${spike_info}_SpikeAlign.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
    samtools sort -@ 25 -o ${align_spike_dir}/${PAIR}_${spike_info}.bam
    samtools index -@ 25 ${align_spike_dir}/${PAIR}_${spike_info}.bam
fi
done

wait

#step 3.3
####### deduplicating with UMIs for experimental genome
echo -e "\n***************************\nRemoving duplicates of experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}

for FILE in ${align_exp_dir}/*.bam
do
    if [ ! -s "${exp_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam)" ]
    then
        (umi_tools dedup \
        -I "$FILE" \
        --umi-separator=":UMI_" \
        --paired \
        -S "${exp_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam)" \
        )> "${rmdup_exp_log}/$(basename ${FILE%.bam}_rmdup.log)" &&
        samtools index "${exp_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam)" 
        samtools flagstat ${exp_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam) > ${rmdup_exp_log}/$(basename ${FILE%.bam}_rmdup.stat)
        echo "rmdupQ10 flagstat has done;file has generated in ${rmdup_exp_log}/$(basename ${FILE%.bam}_rmdup.stat)"
    fi
done

wait


#step 3.4
####### deduplicating with UMIs for spike-in genome
echo -e "\n***************************\nRemoving duplicates of spike-in genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${spike_bam_rmdup}
mkdir -p ${rmdup_spike_log}

for FILE in ${align_spike_dir}/*.bam
do
    if [ ! -s "${spike_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam)" ]
    then
        (umi_tools dedup \
        -I "$FILE" \
        --umi-separator=":UMI_" \
        --paired \
        -S "${spike_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam)" \
        )> "${rmdup_spike_log}/$(basename ${FILE%.bam}_rmdup.log)" &&
        samtools index "${spike_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam)" 
        samtools flagstat "${spike_bam_rmdup}/$(basename ${FILE%.bam}_rmdup.bam)" > ${rmdup_spike_log}/$(basename ${FILE%.bam}_rmdup.stat)
        echo "rmdupQ10 flagstat has done;file has generated in ${rmdup_spike_log}/$(basename ${FILE%.bam}_rmdup.stat)"
    fi
done


wait


#step 3.5
### calculate normalization factors ###
echo -e "\n***************************\nCalculating normalization factors at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

exp_path=${rmdup_exp_log}
spike_path=${rmdup_spike_log}
align_path=${alignexp_log_dir}
scalefactor_file=${logs_dir}/scalefactor_${datainfo}.txt

if [ -s ${scalefactor_file} ] 
then
    rm  ${scalefactor_file}
fi

touch ${scalefactor_file}
ls $align_path|while read file;
do
    sample=${file/_align.log/}
    ALLREADS=$(cat ${align_path}/$file|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    exp_READS=$(cat ${align_path}/$file| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 4,5,8,13,14 |
 sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2+$3)+$4+$5}')
    exp_mapping_RATIO=$(cat ${align_path}/$file|grep "overall alignment rate"|cut -d "%" -f 1)
    exp_qc_READS=$(cat ${exp_path}/${sample}_${exp_info}_rmdup.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
    exp_qc_RATIO=$(echo "${exp_qc_READS}/${ALLREADS}"|bc -l)
    spike_qc_READS=$(cat ${spike_path}/${sample}_${spike_info}_rmdup.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
    QC_reads=$(echo "${exp_qc_READS}+${spike_qc_READS}"|bc )
    spike_qc_RATIO_intotal=$(echo "${spike_qc_READS}/${ALLREADS}"|bc -l)
    spike_qc_RATIO_inqc=$(echo "${spike_qc_READS}/${QC_reads}"|bc -l)
    SCALEFACTOR=$(echo "1000000/${spike_qc_READS}"|bc -l)
    echo -e $sample"\t"$ALLREADS"\t"$exp_READS"\t"$exp_mapping_RATIO"\t"$exp_qc_READS"\t"$exp_qc_RATIO"\t"$spike_qc_READS"\t"$spike_qc_RATIO_intotal"\t"$spike_qc_RATIO_inqc"\t"$SCALEFACTOR >> ${scalefactor_file}
done

wait


#step 4.1
### Making CPM-normalized bigWig files with full-length reads adjusted by spike-in / CPM equal to RPM###
echo -e "\n***************************\nMaking CPM-normalized bigWig files with full-length reads at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${bw_fulllength_dir}

cat  ${scalefactor_file} | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_${exp_info}_rmdup.bam
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_rmdup_fwd.bw" ]
    then
        echo "Generating file: ${bw_fulllength_dir}/${sample}_fulllength_rmdup_fwd.bw "
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_rmdup_fwd.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --samFlagInclude 82
    fi
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_rmdup_rev.bw" ]
    then
        echo "Generating file: ${bw_fulllength_dir}/${sample}_fulllength_rmdup_rev.bw "
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_rmdup_rev.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --samFlagInclude 98
    fi
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_rmdup_rev_minus.bw" ]
    then
        echo "Generating file: ${bw_fulllength_dir}/${sample}_fulllength_rmdup_rev_minus.bw "
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_rmdup_rev_minus.bw \
        --binSize 1 \
        --scaleFactor  -"${scalefactor}" \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --samFlagInclude 98
    fi
done

wait

#step 4.2
### Making CPM-normalized bigWig files with single base adjusted by spike-in for forward and reverse strands###
echo -e "\n***************************\nMaking CPM-normalized bigWig files with single base at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${bw_singlebase_dir}

cat  ${scalefactor_file} | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_${exp_info}_rmdup.bam
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_rmdup_fwd.bw" ]
    then
        echo "Generating file: ${bw_singlebase_dir}/${sample}_singlebase_rmdup_fwd.bw "
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_rmdup_fwd.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 82
    fi
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_rmdup_rev.bw" ]
    then
        echo "Generating file: ${bw_singlebase_dir}/${sample}_singlebase_rmdup_rev.bw "
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_rmdup_rev.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 98
    fi
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_rmdup_rev_minus.bw" ]
    then
        echo "Generating file: ${bw_singlebase_dir}/${sample}_singlebase_rmdup_rev_minus.bw "
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_rmdup_rev_minus.bw \
        --binSize 1 \
        --scaleFactor  -"${scalefactor}" \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 98
    fi
done





