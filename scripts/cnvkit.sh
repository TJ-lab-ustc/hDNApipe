#!/bin/bash
set -e
set -x

# function: call germline or somatic CNV by CNVKIT
# usuage: script_dir/cnvkit.sh $dect_mode $sample_table $out_dir $threads $bin_size $region $seq_mothod


# 1. Import variables 
dect_mode=$1
sample_table=$2
out_dir=$3
threads=$4
bin_size=$5
region=$6
seq_method=$7
dnapipe_dir=$(dirname $(dirname $0))
config="$dnapipe_dir/dnapipe.config"
source $config

eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py3

if [ "$seq_method" != "wgs" ]; then seq_method="amplicon"; fi
if [ -n "$region" ]; then targets="-t ${region}"; fi

# 2. Run CNVKIT
mkdir -p $out_dir/cnvkit/else
mkdir -p $out_dir/cnvkit/plots
# 2-1. somatic 
if [ $dect_mode == "somatic" ]; then
    # get tumor and normal bam
    while IFS=',' read -r name input1 input2 condition sex; do
        if [ $condition == "tumor" ]; then 
            tumor_bam=$out_dir/bam/$name.prcsd.bam
            tumor_name=$name
            if [ -n "$sex" ]; then 
                sample_sex_para="--sample-sex $sex"
            fi
        else
            control_bam=$out_dir/bam/$name.prcsd.bam
            control_name=$name
            if [ "$sex" == "male" ]; then 
                ref_sex_para="-y"
            fi
        fi
    done < $sample_table
    cnvkit.py batch \
        $tumor_bam \
        --normal $control_bam \
        -m $seq_method \
        --fasta $ref \
        $targets \
        --target-avg-size $bin_size \
        --annotate $cnvkit_annotate \
        --drop-low-coverage \
        $ref_sex_para \
        -p $threads \
        --output-reference $out_dir/cnvkit/else/${control_name}.reference.cnn \
        -d $out_dir/cnvkit/else
    cnvkit.py export vcf \
        $out_dir/cnvkit/else/${tumor_name}.prcsd.call.cns \
        -i $tumor_name \
        -o $out_dir/cnvkit/cnv_cnvkit.vcf \
        $sample_sex_para \
        $ref_sex_para
    cnvkit.py scatter \
        -i $tumor_name \
        -o $out_dir/cnvkit/plots/${tumor_name}_scatter.pdf \
        -s $out_dir/cnvkit/else/${tumor_name}.prcsd.cns \
        $out_dir/cnvkit/else/${tumor_name}.prcsd.cnr
    cnvkit.py diagram \
        -o $out_dir/cnvkit/plots/${tumor_name}_diagram.pdf \
        -s $out_dir/cnvkit/else/${tumor_name}.prcsd.cns \
        $out_dir/cnvkit/else/${tumor_name}.prcsd.cnr \
        $sample_sex_para \
        $ref_sex_para \
        --no-gene-labels \
        --title $tumor_name
    cnvkit.py diagram \
        -o $out_dir/cnvkit/plots/${tumor_name}_diagram_label.pdf \
        -s $out_dir/cnvkit/else/${tumor_name}.prcsd.cns \
        $out_dir/cnvkit/else/${tumor_name}.prcsd.cnr \
        $sample_sex_para \
        $ref_sex_para \
        --title $tumor_name
        
# 2-2. germline
else
    mkdir -p $out_dir/cnvkit/vcfs
    mkdir -p $out_dir/cnvkit/call_cns
    mkdir -p $out_dir/cnvkit/plots
    while IFS=',' read -r name input1 input2 condition sex; do
        bam=$out_dir/bam/$name.prcsd.bam
        if [ -n "$sex" ]; then sample_sex_para="--sample-sex $sex"; fi
        if [ "$sex" == "male" ]; then ref_sex_para="-y"; fi
        cnvkit.py batch \
            $bam \
            -n \
            -m $seq_method \
            --fasta $ref \
            $targets \
            --target-avg-size $bin_size \
            --annotate $cnvkit_annotate \
            $ref_sex_para \
            -p $threads \
            --output-reference $out_dir/cnvkit/else/${name}.reference.cnn \
            -d $out_dir/cnvkit/else
        cnvkit.py export vcf \
            -o $out_dir/cnvkit/vcfs/${name}_cnvkit.vcf \
            $out_dir/cnvkit/else/$name.prcsd.call.cns \
            $ref_sex_para \
            $sex_para \
            --sample-id $name
        cnvkit.py scatter \
            -i $name \
            -o $out_dir/cnvkit/plots/${name}_scatter.pdf \
            -s $out_dir/cnvkit/else/${name}.prcsd.cns \
            $out_dir/cnvkit/else/${name}.prcsd.cnr
    done < $sample_table
    cnvkit.py heatmap \
        -o $out_dir/cnvkit/plots/samples_heatmap.pdf \
        $out_dir/cnvkit/else/*.prcsd.cns \
        -d
    mv $out_dir/cnvkit/else/*.prcsd.call.cns $out_dir/cnvkit/call_cns/
fi
