#!/usr/bin/bash

# function
function CheckBam()
{
    bam=$1
    if [ -s $bam ]
    then
        samtools quickcheck $bam
         if [ $? -eq 0 ]
            then
                touch ${bam}.done
            else
                echo "${bam} is truncked, ${bam} is deleted"
                rm $bam
            fi
        fi
            
}


# usage
if [ $# -ne 3 ]; then
    echo -e $0: usage: bash $0 "$1: Path of directory; $2: datainfo(eg: RNA_220701) and $3:reference genome"
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
trimmedFastq_log_dir=${file_path}/logs/trimmeddata
trimmedFastq_fastqc_dir=${file_path}/02.2_trimmeddata_fastqc

align_exp_dir=${file_path}/03_bam
alignexp_log_dir=${file_path}/logs/align_exp
align_spike_dir=${file_path}/03_spikebam
alignspike_log_dir=${file_path}/logs/align_spike
exp_bam_rmdup=${file_path}/03_bam_rmdup
rmdup_exp_log=${file_path}/logs/rmdup_state_exp
spike_bam_rmdup=${file_path}/03_spikebam_rmdup
rmdup_spike_log=${file_path}/logs/rmdup_state_spike

bw_track_dir=${file_path}/04_bigwig_track

echo -e "*RNA-seq parameters:
          Run_path: ${file_path}
          data_path: ${input_path}
          reference_genome: ${reference_info}"

#reference genome
GENOME_human="/share/home/Blueberry/reference/index/star/hg19_star_index"
GENOME_mouse="/share/home/Blueberry/reference/index/star/mm10_star_index"

if [[ ${reference_info} == "hg19" ]]
then
    GENOME_EXP=${GENOME_human}
    GENOME_SPIKE=${GENOME_mouse}
    exp_info="hg19"
    spike_info="mm10"
else
    GENOME_EXP=${GENOME_mouse}
    GENOME_SPIKE=${GENOME_human}
    exp_info="mm10"
    spike_info="hg19"
fi

echo -e "\n***************************\nTT-seq processing at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
echo -e "experimental genome is: ${GENOME_EXP} \nspike-in genome is: ${GENOME_SPIKE}"


#================== rawdir ==================#
### rename if exists sampleinfo.txt
#step 1.1 change file name#####
echo -e "\n***************************\nRenaming files at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${raw_dir}

cat $sampleinfo| while read id;
do
    arr=($id)
    sample1=${arr[0]}
    sample2=${arr[1]}
    if [ ! -s ${raw_dir}/${sample2}_R1.fastq.gz ]
    then
        fq1=$(ls ${input_path}/*/*1.f*q.gz|grep "$sample1")
        fq2=$(ls ${input_path}/*/*2.f*q.gz|grep "$sample1")
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
multiqc ${fastqc_dir}/ -n rawdata_multiqc_${datainfo} -o ${fastqc_dir}/


#================== trim ==================#
#step 2.1
### Trimming adapters  (trim_galore)

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
        trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 --paired -o ${trimmedFastq_dir} $fq1 $fq2 \
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
#####qc for trimmed data####
echo -e "\n***************************\nqc for trimmed data at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${trimmedFastq_fastqc_dir}

nohup_number=0
for FILE in `ls ${trimmedFastq_dir}/*fq.gz`
do
    sample=$(basename ${FILE/.fq.gz/})
    if [ ! -s ${trimmedFastq_fastqc_dir}/"$(basename ${FILE/.fq.gz/_fastqc.html})" ]
    then
        echo "Generating file: ${trimmedFastq_fastqc_dir}/"$(basename ${FILE/.fq.gz/_fastqc.html})"; "
        fastqc $FILE -t 4 -o ${trimmedFastq_fastqc_dir} >  ${trimmedFastq_fastqc_dir}/${sample}_fastqc.log 2>&1 &
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



#================== align ==================#
#step 3.1
###Aligning to experimental genome#####
echo -e "\n***************************\nAligning to experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
mkdir -p  ${align_exp_dir}
mkdir -p  ${alignexp_log_dir}


for PAIR in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | uniq )
do
    CheckBam ${align_exp_dir}/${PAIR}_${exp_info}.Aligned.sortedByCoord.out.bam
    if [ ! -s ${align_exp_dir}/${PAIR}_${exp_info}.Aligned.sortedByCoord.out.bam ]
    then
        STAR --runThreadN 27 --genomeDir ${GENOME_EXP} \
        --readFilesCommand zcat --readFilesIn ${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz ${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterType BySJout --outFilterMultimapNmax 1 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFileNamePrefix ${align_exp_dir}/${PAIR}_${exp_info}.
    fi
done


#step 3.2
###Aligning to spike-in genome to get normalization factors#####
echo -e "\n***************************\nAligning to spike-in genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${align_spike_dir}
mkdir -p ${alignspike_log_dir}


for PAIR in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | uniq )
do
    CheckBam ${align_spike_dir}/${PAIR}_${spike_info}.Aligned.sortedByCoord.out.bam
    if [ ! -s ${align_spike_dir}/${PAIR}_${spike_info}.Aligned.sortedByCoord.out.bam ]
    then
        STAR --runThreadN 27 --genomeDir ${GENOME_SPIKE} \
        --readFilesCommand zcat --readFilesIn ${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz ${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterType BySJout --outFilterMultimapNmax 1 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFileNamePrefix ${align_spike_dir}/${PAIR}_${spike_info}.
    fi
done



#================== rmdup ==================#
#step 3.3 & 3.4
####### deduplicating with UMIs for experimental genome
echo -e "\n***************************\nRemoving duplicates of experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}

mkdir -p ${spike_bam_rmdup}
mkdir -p ${rmdup_spike_log}


for i in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | uniq )
do
    if [ ! -s ${rmdup_exp_log}/${i}_${exp_info}.rmdup.q7.stat ]
    then

    input=`ls ${align_exp_dir}/${i}_${exp_info}*.sortedByCoord.out.bam`
    rmdup=${exp_bam_rmdup}/${i}_${exp_info}.rmdup.bam
    metrics=${exp_bam_rmdup}/${i}_${exp_info}.metrics
    picard MarkDuplicates REMOVE_DUPLICATES=True \
        INPUT=$input OUTPUT=${rmdup} METRICS_FILE=${metrics} \
        2>${rmdup_exp_log}/${i}_${exp_info}.dup.log;
    rmdupQ7=${exp_bam_rmdup}/${i}_${exp_info}.rmdup.q7.bam
    samtools view -b -q 7 -f 2 ${rmdup} -o ${rmdupQ7}
    samtools index -@ 26 ${rmdupQ7}
    samtools flagstat -@ 25 ${rmdupQ7} > ${rmdup_exp_log}/${i}_${exp_info}.rmdup.q7.stat
    echo "Exp: rmdupQ7 flagstat has done;file has generated in ${rmdup_exp_log}/${i}_${exp_info}.rmdup.q7.stat"
    else
    echo "${rmdup_exp_log}/${i}_${exp_info}.rmdup.q7.stat exists, continue..."
    fi

    if [ ! -s ${rmdup_spike_log}/${i}_${spike_info}.rmdup.q7.stat ]
    then
    input=`ls ${align_spike_dir}/${i}_${spike_info}*.sortedByCoord.out.bam`
    rmdup=${spike_bam_rmdup}/${i}_${spike_info}.rmdup.bam
    metrics=${spike_bam_rmdup}/${i}_${spike_info}.metrics
    picard MarkDuplicates REMOVE_DUPLICATES=True \
        INPUT=$input OUTPUT=${rmdup} \
        METRICS_FILE=${metrics} 2>${rmdup_spike_log}/${i}_${spike_info}.dup.log;
    rmdupQ7=${spike_bam_rmdup}/${i}_${spike_info}.rmdup.q7.bam
    samtools view -b -q 7 -f 2 ${rmdup} -o ${rmdupQ7}
    samtools index -@ 25 ${rmdupQ7}
    samtools flagstat -@ 25 ${rmdupQ7} > ${rmdup_spike_log}/${i}_${spike_info}.rmdup.q7.stat
    echo "Spike: rmdupQ7 flagstat has done;file has generated in ${rmdup_spike_log}/${i}_${spike_info}.rmdup.q7.stat"
    else
    echo "${rmdup_spike_log}/${i}_${spike_info}.rmdup.q7.stat exists, continue..."
    fi
done


#================== spikein ==================#
#step 3.5
### calculate normalization factors ###
echo -e "\n***************************\nCalculating normalization factors at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

exp_path=${rmdup_exp_log}
spike_path=${rmdup_spike_log}

scalefactor_file=${logs_dir}/scalefactor_${datainfo}.txt

echo -e "sample\tALLREADS\tExp_READS\tExp_RATIO\tExp_rmdupQ7_READS\tExp_rmdupQ7_RATIO\tSpikeIn_READS\tSpikeIn_RATIO\tSpikeIn_rmdupQ7_READS\tSpikeIn_rmdupQ7_RATIO\tSCALEFACTOR" >${scalefactor_file}

for file in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | uniq );
do
sample=$(basename $file)  #${file/_align.log/}
ALLREADS=`grep "Number of input reads" ${align_exp_dir}/${sample}_${exp_info}*Log.final.out | awk '{print $NF*2}'`
Exp_READS=`grep -e "Number.* mapped" -e "mapped.*number" ${align_exp_dir}/${sample}_${exp_info}*Log.final.out | cut -f 2 |  xargs | awk '{print ($1+$2+$3)*2}'`
Exp_RATIO=`grep " mapped.*%" ${align_exp_dir}/${sample}_${exp_info}*Log.final.out | cut -f 2 | xargs | awk '{print $1+$2+$3}'`
Exp_rmdupQ7_READS=$(grep "total (QC-passed reads" ${exp_path}/${sample}_${exp_info}.rmdup.q7.stat |awk '{print $1}')
Exp_rmdupQ7_RATIO=$(echo "${Exp_rmdupQ7_READS}/${ALLREADS}"|bc -l)
SpikeIn_READS=`grep -e "Number.* mapped" -e "mapped.*number" ${align_spike_dir}/${sample}_${spike_info}*Log.final.out | cut -f 2 |  xargs | awk '{print ($1+$2+$3)*2}'`
SpikeIn_RATIO=`grep " mapped.*%" ${align_spike_dir}/${sample}_${spike_info}*Log.final.out | cut -f 2 | xargs | awk '{print $1+$2+$3}'`
SpikeIn_rmdupQ7_READS=$(grep "total (QC-passed reads" ${spike_path}/${sample}_${spike_info}.rmdup.q7.stat |cut -d " " -f 1)
SpikeIn_rmdupQ7_RATIO=$(echo "${SpikeIn_rmdupQ7_READS}/${ALLREADS}"|bc -l)
SCALEFACTOR=$(echo "1000000/$SpikeIn_rmdupQ7_READS"|bc -l)
echo -e ${sample}"\t"$ALLREADS"\t"$Exp_READS"\t"$Exp_RATIO"\t"$Exp_rmdupQ7_READS"\t"$Exp_rmdupQ7_RATIO"\t"$SpikeIn_READS"\t"$SpikeIn_RATIO"\t"$SpikeIn_rmdupQ7_READS"\t"$SpikeIn_rmdupQ7_RATIO"\t"$SCALEFACTOR >> ${scalefactor_file}
done

wait


#================== bamcoverage(bw) ==================#
echo -e "\n
***************************
track of R1andR2 begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"

mkdir -p ${bw_track_dir}

sed 1d ${scalefactor_file} | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[10]}
bam_file=${exp_bam_rmdup}/${sample}_${exp_info}.rmdup.q7.bam
echo "=== sample: $sample ==="
    if [ ! -s "${bw_track_dir}/${sample}_plus.bw" ]
    then
        echo "Generating file: ${bw_track_dir}/${sample}_plus.bw..."
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_track_dir}/${sample}_plus.bw \
        --binSize 1 \
        --scaleFactor $scalefactor \
        --numberOfProcessors 25 \
        --normalizeUsing None \
        --filterRNAstrand forward
    else
        echo "${bw_track_dir}/${sample}_plus.bw exists, continue..."
    fi

    if [ ! -s "${bw_track_dir}/${sample}_rev.bw" ]
    then
    echo "Generating file: ${bw_track_dir}/${sample}_rev.bw"
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_track_dir}/${sample}_rev.bw \
        --binSize 1 \
        --scaleFactor $scalefactor \
        --numberOfProcessors 25 \
        --normalizeUsing None \
        --filterRNAstrand reverse
    else
    echo "${bw_track_dir}/${sample}_rev.bw exists. continue..."
    fi

    if [ ! -s "${bw_track_dir}/${sample}_minus.bw" ]
    then
    echo "Generating file: ${bw_track_dir}/${sample}_minus.bw..."
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_track_dir}/${sample}_minus.bw \
        --binSize 1 \
        --scaleFactor -"${scalefactor}" \
        --numberOfProcessors 25 \
        --normalizeUsing None \
        --filterRNAstrand reverse
    else
    echo "${bw_track_dir}/${sample}_minus.bw exists, continue"
    fi
done
