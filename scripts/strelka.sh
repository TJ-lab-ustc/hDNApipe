#!/bin/bash
set -e
set -o pipefail

# function: call germline or somatic SNP/INDEL by strelka2
# usuage: DIR/strelka.sh $dect_mode $seq_method $sample_table $out_dir $threads $region
# example: ./strelka.sh germline-cohort wgs /folder/sample.txt ./ 20 ''

# 1. Import variables 
dect_mode=$1
seq_method=$2
sample_table=$3
out_dir=$(readlink -f $4)
threads=$5
region=$6

if [ $seq_method == "wes" ]; then 
    wes_para="--exome"
fi
if [ -n "$region" ]; then
    region_para="--callRegions ${region}.gz"
fi

# 2. Import config
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

# 3. Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py2

# 4. INDEX USER BED REGION
if [ -f "$region" ] && [ ! -f ${region}.gz.tbi ]; then
    conda run -n $dnapipe_env_py3 bgzip -fk $region
    conda run -n $dnapipe_env_py3 tabix -f -p bed ${region}.gz
fi

# 5. Run strelka2
strelka_dir=$out_dir/strelka
mkdir -p $strelka_dir
BCFTOOLS="conda run -n $dnapipe_env_py3 bcftools"

# 5-1 somatic
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
    configureStrelkaSomaticWorkflow.py \
        --tumorBam $tumor_bam \
        --normalBam $control_bam \
        --referenceFasta $ref \
        --runDir $strelka_dir \
        $region_para \
        $wes_para
    # run step 2
    $strelka_dir/runWorkflow.py \
        -m local \
        -j $threads
    # extract
    snp_vcf=$strelka_dir/results/variants/somatic.snvs.vcf.gz
    indel_vcf=$strelka_dir/results/variants/somatic.indels.vcf.gz
    $BCFTOOLS view -f PASS $snp_vcf -o $strelka_dir/snv_strelka.vcf -O v
    $BCFTOOLS view -f PASS $indel_vcf -o $strelka_dir/indel_strelka.vcf -O v

# 5-2 germline
else
    # get bam list
    bam_list=""
    while IFS=',' read -r name input1 input2 condition sex; do
        bam=$out_dir/bam/$name.prcsd.bam
        bam_list+="--bam $bam "
    done < $sample_table
    # get scripts
    configureStrelkaGermlineWorkflow.py \
        $bam_list \
        --referenceFasta $ref \
        --runDir $strelka_dir \
        $region_para \
        $wes_para
    # run step 2
    $strelka_dir/runWorkflow.py \
        -m local \
        -j $threads
    # output file
    strelka_vcf=$strelka_dir/results/variants/variants.vcf.gz
    # extract
    $BCFTOOLS view \
        -f PASS \
        -v snps \
        -o $strelka_dir/snv_strelka.vcf \
        $strelka_vcf
    $BCFTOOLS view \
        -f PASS \
        -v indels \
        -o $strelka_dir/indel_strelka.vcf \
        $strelka_vcf
fi

