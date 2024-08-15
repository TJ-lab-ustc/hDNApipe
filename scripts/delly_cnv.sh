#!/bin/bash
set -e

# function: call germline or somatic CNV by DELLY CNV (only for WGS data)
# usuage: script_dir/delly_cnv.sh $dect_mode $sample_table $out_dir $bin_size
# note: only work on WGS data; no any region parameter

# 1. Import variables 
dect_mode=$1
sample_table=$2
out_dir=$3
bin_size=$4
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py3


# 2. Run DELLY
mkdir -p $out_dir/delly_cnv/else
mkdir -p $out_dir/delly_cnv/figure
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
    echo -e "${tumor_name}\ttumor\n${control_name}\tcontrol\n" > $out_dir/delly_cnv/else/delly_sample.tsv
    # run
    # delly cnv \
    #     -u \
    #     -z $bin_size \
    #     -o $out_dir/delly_cnv/else/tumor.bcf \
    #     -c $out_dir/delly_cnv/else/tumor.cov.gz \
    #     -g $ref \
    #     -m $map \
    #     $tumor_bam
    # delly cnv \
    #     -u \
    #     -v $out_dir/delly_cnv/else/tumor.bcf \
    #     -o $out_dir/delly_cnv/else/control.bcf \
    #     -g $ref \
    #     -m $map \
    #     $control_bam
    bcftools merge \
        -m id \
        -O b \
        -o $out_dir/delly_cnv/else/tumor_control.bcf \
        $out_dir/delly_cnv/else/tumor.bcf \
        $out_dir/delly_cnv/else/control.bcf
    bcftools index -f $out_dir/delly_cnv/else/tumor_control.bcf
    delly classify \
        -p \
        -f somatic \
        -o $out_dir/delly_cnv/else/somatic.bcf \
        -s $out_dir/delly_cnv/else/delly_sample.tsv \
        $out_dir/delly_cnv/else/tumor_control.bcf
    # bcf to vcf
    bcftools view $out_dir/delly_cnv/else/somatic.bcf > $out_dir/delly_cnv/cnv_delly.vcf
    # plot
    bcftools query \
        -s $tumor_name \
        -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" \
        $out_dir/delly_cnv/else/somatic.bcf \
        > $out_dir/delly_cnv/else/segmentation.bed
    cd $out_dir/delly_cnv/figure
    /usr/bin/Rscript $DELLY_plot_dir/rd.R \
        $out_dir/delly_cnv/else/tumor.cov.gz \
        $out_dir/delly_cnv/else/segmentation.bed

# 2-2. germline
else
    bcf=''
    geno=''
    # Call CNVs for each sample
    while IFS=',' read -r name input1 input2 condition sex; do
        bam=$out_dir/bam/$name.prcsd.bam
        delly cnv \
            -o $out_dir/delly_cnv/else/$name.bcf \
            -g $ref \
            -m $map \
            $bam
        bcfs+="$out_dir/delly_cnv/else/$name.bcf "
    done < $sample_table
    # Merge CNVs into a unified site list
    delly merge \
        -e \
        -p \
        -o $out_dir/delly_cnv/else/sites.bcf \
        $bcfs
    # Genotype CNVs for each sample
    while IFS=',' read -r name input1 input2 condition sex; do
        bam=$out_dir/bam/$name.prcsd.bam
        delly cnv \
            -u \
            -v $out_dir/delly_cnv/else/sites.bcf \
            -g $ref \
            -m $map \
            -o $out_dir/delly_cnv/else/$name.geno.bcf \
            $bam
        geno+="$out_dir/delly_cnv/else/$name.geno.bcf "
        delly_bcf="$out_dir/delly_cnv/else/$name.geno.bcf"
    done < $sample_table
    # Merge genotypes using bcftools
    sample_num=$( cat $sample_table |wc -l )
    if [ $sample_num -gt 1 ]; then
        bcftools merge \
            -m id \
            -O b \
            -o $out_dir/delly_cnv/else/merged.bcf \
            $geno
        delly_bcf=$out_dir/delly_cnv/else/merged.bcf
    fi
    # Filter for germline CNVs
    bcftools index -f $delly_bcf
    delly classify \
        -f germline \
        -o $out_dir/delly_cnv/else/filtered.bcf \
        $delly_bcf
    # Bcf to vcf
    bcftools view $out_dir/delly_cnv/else/filtered.bcf \
        > $out_dir/delly_cnv/cnv_delly.vcf
    # plot
    if [ $sample_num -ge 3 ]; then
        bcftools query \
            -f "%ID[\t%RDCN]\n" \
            $out_dir/delly_cnv/else/filtered.bcf \
            > $out_dir/delly_cnv/else/plot.tsv
        cd $out_dir/delly_cnv/figure
        /usr/bin/Rscript $DELLY_plot_dir/cnv.R $out_dir/delly_cnv/else/plot.tsv
    fi
fi
