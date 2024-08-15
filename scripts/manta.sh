#!/bin/bash
set -e

# function: call germline/somatic short variants(SNP/INDEL) by strelka2
# usuage: DIR/strelka.sh $dect_mode $seq_method $sample_table $out_dir $threads $region

# 1. Import variables 
dect_mode=$1
seq_method=$2
sample_table=$3
out_dir=$4
threads=$5
region=$6
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

# 2. Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py2

# 3. INDEX USER BED REGION
if [ -f "$region" ] && [ ! -f ${region}.gz.tbi ]; then
    conda run -n $dnapipe_env_py3 bgzip -fk $region
    conda run -n $dnapipe_env_py3 tabix -f -p bed ${region}.gz
fi

# 4. Run manta step1
if [ $seq_method == "wes" ]; then 
    wes_para="--exome"
fi
if [ -n "$region" ]; then
    region_para="--callRegions ${region}.gz"
fi
manta_dir=$out_dir/manta/
mkdir -p $manta_dir
BCFTOOLS="conda run -n $dnapipe_env_py3 bcftools"

# 4-1 somatic
if [ $dect_mode == "somatic" ]; then
    # get tumor and normal bam
    while IFS=',' read -r name input1 input2 condition sex; do
        if [ $condition == "tumor" ]; then 
            tumor_bam=$out_dir/bam/$name.prcsd.bam
        else
            control_bam=$out_dir/bam/$name.prcsd.bam
        fi
    done < $sample_table
    # get scripts
    configManta.py \
        --tumorBam $tumor_bam \
        --normalBam $control_bam \
        --referenceFasta $ref \
        --runDir $manta_dir \
        $region_para \
        $wes_para
    manta_vcf="$manta_dir/results/variants/somaticSV.vcf.gz"
    
# 4-2 germline
else
    # get bam list
    while IFS=',' read -r name input1 input2 condition sex; do
        bam=$out_dir/bam/$name.prcsd.bam
        bam_list+="--bam $bam "
    done < $sample_table
    # get scripts
    configManta.py \
        $bam_list \
        --referenceFasta $ref \
        --runDir $manta_dir \
        $region_para \
        $wes_para
    manta_vcf="$manta_dir/results/variants/diploidSV.vcf.gz"
fi

# 5. Run manta step2
$manta_dir/runWorkflow.py -m local -j $threads

# 6. Filter PASS
$BCFTOOLS view -f PASS $manta_vcf -o $manta_dir/sv_manta.vcf -O v
