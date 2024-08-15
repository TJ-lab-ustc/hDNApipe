#!/bin/bash
set -e
set -o pipefail

# function: alignment and/or pre-process of input fastq or bam file
# usuage: $DIR/align_pre.sh $file_type $input1 $input2 $threads $memory $outdir $RG $rm_dup $sup_clip $map_qual "$force_rg"

# 1. bwa alignment
# 2. sort by SAMTOOLS 
# 3. mark duplicates by SAMBAMBA
# for bam: check if sorted or marked duplicates; if not, then run pre-process

# Output files:
# 1. Sort: $out_dir/bam/$sample_name.sorted.bam (if input not sorted)
# 2. Mark Dup: $out_dir/bam/$sample_name.prcsd.bam (soft link if already marked)


# 1. GET VARIABLES
file_type=$1
input1=$2
input2=$3
threads=$4
memory=$5
out_dir=$(readlink -f $6)
RG=$7
rm_dup=$8
sup_clip=$9
map_qual=${10}
force_rg=${11}

if [ "$rm_dup" == "True" ]; then dup_para="--remove-duplicates"; fi
if [ "$sup_clip" == "True" ]; then clip_para="-Y"; fi

# 2. Import config
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

# 3. Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py3

# 4. RUN
bam_dir=$out_dir/bam
temp_dir=$out_dir/temp
mkdir -p $bam_dir
mkdir -p $temp_dir
sample=${RG#*SM:}
# 4-1. ALIGN and SORT
if [ $file_type == "fastq" ];then
    echo "============= $sample align and sort ============="
    bwa mem -R $RG -T $map_qual $clip_para -t $threads $ref $input1 $input2 | \
        samtools sort -@ $threads -m $memory -o $bam_dir/$sample.sorted.bam -
    sort_bam=$bam_dir/$sample.sorted.bam
else
    echo "============= $sample bam preprocessing ============="
    # check if RG information is contained
    samtools view -H $input1 | grep -q "^@RG" || RG_flag=False
    if [ "$RG_flag" == "False" ] || [ "$force_rg" == "True" ]; then
        echo "Adding read group information..."
        samtools addreplacerg -r $RG -@ $threads -m overwrite_all -w -o $bam_dir/${sample}_addRG.bam $input1
        input1="$bam_dir/${sample}_addRG.bam"
    fi
    # sort
    samtools view -H $input1 | grep -q 'SO:coordinate' || sort_flag=False
    if [ "$sort_flag" == "False" ]; then
        echo "Sorting now ..."
        samtools sort -@ $threads -m $memory -o $bam_dir/$sample.sorted.bam $input1
        sort_bam=$bam_dir/$sample.sorted.bam
    else
    # check if marked duplicates
        echo "The bam file is already sorted. Skipping sortment."
        samtools view -H $input1 | grep 'PG' | grep -qE 'CL:MarkDuplicates|CL:markdup' || mark_flag=False
        if [ "$mark_flag" == "False" ]; then
            sort_bam=$input1
        elif [ $input1 == "$bam_dir/$sample.prcsd.bam" ]; then
            echo "The duplicates of the bam file is already marked. Skipping."
        else
            echo "The duplicates of the bam file is already marked. Skipping."
            ln -sf $input1 $bam_dir/$sample.prcsd.bam
        fi
    fi
fi

# 4-2. MARK DUPLICATES
if [ -n "$sort_bam" ]; then
    echo "============= $sample mark duplicates ============="
    ulimit -u 8000
    sambamba markdup -t $threads $dup_para --tmpdir $temp_dir $sort_bam $bam_dir/$sample.prcsd.bam
fi

# 4-3. INDEX BAM
if [ ! -f "$bam_dir/$sample.prcsd.bam.bai" ]; then
    echo "============= $sample bam file index ============="
    sambamba index -q -t $threads $bam_dir/$sample.prcsd.bam $bam_dir/$sample.prcsd.bam.bai
fi

echo "============= $sample pre-process done ============="

