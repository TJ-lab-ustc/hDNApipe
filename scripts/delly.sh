#!/bin/bash
set -e

# function: call germline or somatic SV by delly (only for WGS data)
# usuage: script_dir/delly.sh $dect_mode $sample_table $out_dir $region
# note: only work on WGS data; default exclude specific regions such as telomere; 
#       no other region parameters.


# 1. Import variables
dect_mode=$1
sample_table=$2
out_dir=$3
region=$4
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

region_para="-x $excl_region"

eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py3

# 2. Run DELLY
mkdir -p $out_dir/delly/else
# 2-1. somatic 
if [ $dect_mode == "somatic" ]; then
    # get tumor and normal bam
    while IFS=',' read -r name input1 input2 condition sex; do
        if [ $condition == "tumor" ]; then 
            tumor_bam=$out_dir/bam/$name.prcsd.bam
            tumor_name=$name
        else
            control_bam=$out_dir/bam/$name.prcsd.bam
            control_name=$name
        fi
    done < $sample_table
    # get delly tsv
    echo -e "${tumor_name}\ttumor\n${control_name}\tcontrol\n" > $out_dir/delly/else/delly_sample.tsv
    # run
    delly call \
        $region_para \
        -o $out_dir/delly/else/orig.bcf \
        -g $ref \
        $tumor_bam \
        $control_bam
    delly filter \
        -p \
        -f somatic \
        -o $out_dir/delly/else/pre.bcf \
        -s $out_dir/delly/else/delly_sample.tsv \
        $out_dir/delly/else/orig.bcf
    bcftools view $out_dir/delly/else/pre.bcf > $out_dir/delly/sv_delly.vcf
# 2-2. germline
else
    bcf=''
    geno=''
    # call bcf
    while IFS=',' read -r name input1 input2 condition sex; do
        bam=$out_dir/bam/$name.prcsd.bam
        delly call \
            -g $ref \
            $region_para \
            -o $out_dir/delly/else/$name.bcf \
            $bam
        bcf=$out_dir/delly/else/$name.bcf
        bcfs+="$out_dir/delly/else/$name.bcf "
    done < $sample_table
    sample_num=$( cat $sample_table |wc -l )
    if [ $sample_num -gt 1 ]; then
        delly merge \
            -p \
            -o $out_dir/delly/else/merge_sites.bcf \
            $bcfs
        bcf=$out_dir/delly/else/merge_sites.bcf
    fi
    # genotype
    while IFS=',' read -r name input1 input2 condition sex; do
        delly call \
            $region_para \
            -g $ref \
            -v $bcf \
            -o $out_dir/delly/else/$name.geno.bcf \
            $out_dir/bam/$name.prcsd.bam
        geno+="$out_dir/delly/else/$name.geno.bcf "
        delly_bcf="$out_dir/delly/else/$name.geno.bcf"
    done < $sample_table
    if [ $sample_num -gt 1 ]; then
        bcftools merge \
            -m id \
            -o $out_dir/delly/else/merge.bcf \
            $geno
        delly_bcf="$out_dir/delly/else/merge.bcf"
    fi
    # delly group variants filter
    if [ $sample_num -ge 20 ]; then 
        delly filter \
            -f germline \
            -o $out_dir/delly/else/germline.bcf \
            $out_dir/delly/else/merge.bcf
        delly_bcf="$out_dir/delly/else/germline.bcf"
    fi
    # bcf to vcf
    bcftools view \
        -e '( GT[*]="mis" || GT[*]="0/0")' \
        $delly_bcf \
        -o $out_dir/delly/sv_delly.vcf
fi
