#!/bin/bash
set -e

# Usage: bash annotate.sh $vcf_file $out_dir $threads [annotation_paramters]
# Function: annotate vcf by VEP and AnnotSV

# 1. Import variables 
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config
vcf=$1
out_dir=$2
threads=$3
mkdir -p $out_dir/annot
stats_para="--no_stats"


eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py3

for para in $@; do
    case $para in
        "exist")
            exist_para="--check_existing"
        ;;
        "symbol")
            symbol_para="--symbol"
        ;;
        "report")
            report_flag="True"
            stats_para=""
        ;;
        "pathogenicity")
            snv_patho_para="--sift b --polyphen b"
            indel_patho_para="--plugin CADD,$cadd_indel_file"
        ;;
        "frequency")
            af_para="--af --af_1kg --af_gnomade --af_gnomadg --max_af"
        ;;
        "clinvar")
            clinvar_para="--custom $clinvar_file,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN"
        ;;
    esac
done

# 2. AnnotSV for SV/CNV
name=$(basename $vcf)
vcf_prefix=${name%.vcf}
if [[ $(echo "$name" | grep 'sv') ]] || [[ $(echo "$name" | grep 'cnv') ]]; then
    AnnotSV \
        -SVinputFile $vcf \
        -genomeBuild GRCh38 \
        -outputFile $out_dir/annot/${vcf_prefix}_annotsv \
        -annotationsDir $annotsv_dir \
        -vcf 1
    out_tsv=$out_dir/annot/${vcf_prefix}_annotsv.tsv
    rm $out_dir/annot/*.log
    if [ "$report_flag" == "True" ]; then
        mkdir -p $out_dir/annot/html
        perl $KNOTANNOTSV \
            --annotSVfile $out_tsv \
            --outDir $out_dir/annot/html \
            -genomeBuild GRCh38 \
            --configFile $config_knot
    fi

# 3. VEP for snp/indel
elif [[ $(echo "$name" | grep 'snv') ]] || [[ $(echo "$name" | grep 'indel') ]]; then
    # SIFT/PolyPhen for SNV, CADD for INDEL
    if [[ $(echo "$name" | grep 'snv') ]]; then
        indel_patho_para=""
    else
        snv_patho_para=""
    fi
    # RUN
    vep \
        --vcf --force_overwrite --cache \
        -i $vcf \
        --fork $threads \
        -o $out_dir/annot/${vcf_prefix}_vep.vcf \
        --dir_plugins ${dir_plugins} \
        --dir_cache ${dir_cache} \
        $exist_para $symbol_para $stats_para $snv_patho_para $af_para $clinvar_para $indel_patho_para
    if [ "$report_flag" == "True" ]; then
        mkdir -p $out_dir/annot/html
        mv $out_dir/annot/*.html $out_dir/annot/html
    fi
fi

