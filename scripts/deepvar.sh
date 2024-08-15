#!/bin/bash
set -e

# function: call germline short variants(SNP/INDEL) by deepvariant
# usuage: DIR/deepvar.sh $seq_method $sample_table $out_dir $threads $region
# example: ./deepvar.sh wgs /folder/sample.txt ./ 20 ''

# 1. Import variables 
seq_method=$1
sample_table=$2
out_dir=$(readlink -f $3)
threads=$4
region=$5

if [ $seq_method == "wgs" ]; then 
    model_para="WGS"
    gvcf_para="DeepVariantWGS"
else
    model_para="WES"
    gvcf_para="DeepVariantWES"
    gvcf_bed="--bed $region"
fi

# 2. Import config
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

# 3. Run deepvariant in main environment
deepvar_dir=$out_dir/deepvar
mkdir -p $deepvar_dir/gvcf
mkdir -p $deepvar_dir/vcfs
mkdir -p $deepvar_dir/html
gvcfs=''
while IFS=',' read -r name input1 input2 condition sex; do
    bam=$out_dir/bam/$name.prcsd.bam
    run_deepvariant \
        --model_type $model_para \
        --ref $ref \
        --reads $bam \
        --regions $region \
        --output_vcf $deepvar_dir/vcfs/$name.output.vcf.gz \
        --output_gvcf $deepvar_dir/gvcf/$name.output.g.vcf.gz \
        --num_shards $threads
    deepvar_vcf=$deepvar_dir/vcfs/$name.output.vcf.gz
    gvcfs+="$deepvar_dir/gvcf/$name.output.g.vcf.gz "
done < $sample_table

# 4. Change conda environment
eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py3

# 5. Combine gvcfs
sample_num=$( less $sample_table |wc -l )
if [ $sample_num -gt 1 ]; then
    glnexus_cli \
        --dir $deepvar_dir/GLnexus.DB \
        --config $gvcf_para \
        --threads $threads \
        $gvcf_bed \
        $gvcfs | \
            bcftools view - | \
            bgzip -c \
            > $deepvar_dir/vcfs/deepvariant.cohort.vcf.gz
    deepvar_vcf=$deepvar_dir/vcfs/deepvariant.cohort.vcf.gz
fi

# 5. Extract
if [ $sample_num -eq 1 ]; then
    pass_para="-f PASS"
fi
bcftools view \
    $pass_para \
    -v snps \
    -o $deepvar_dir/snv_deepvar.vcf \
    $deepvar_vcf
bcftools view \
    $pass_para \
    -v indels \
    -o $deepvar_dir/indel_deepvar.vcf \
    $deepvar_vcf
mv $deepvar_dir/vcfs/*.html $deepvar_dir/html