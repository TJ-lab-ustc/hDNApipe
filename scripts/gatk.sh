#!/bin/bash
set -e

# function: call germline or somatic SNP/INDEL by gatk HaplotypeCaller or Mutect2
# usuage: script_dir/Mutect2.sh $dect_mode $seq_method $sample_table $out_dir $region_bed

# 1. Import variables 
dect_mode=$1
seq_method=$2
sample_table=$3
out_dir=$4
region=$5
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

if [ -n "$region" ]; then
    interval_para="-L $region"
fi

BCFTOOLS="conda run -n $dnapipe_env_py3 bcftools"

# 2. Check dict
ref_prefix=${ref%.fa}
ref_prefix=${ref_prefix%.fa.gz}
ref_prefix=${ref_prefix%.fasta}
ref_prefix=${ref_prefix%.fasta.gz}
ref_prefix=${ref_prefix%.fna}
ref_prefix=${ref_prefix%.fna.gz}
if [ ! -f "$ref_prefix.dict" ]; then
    echo "**********Reference dictionary not created. Creating now...**********"
    $GATK CreateSequenceDictionary -R $ref
    echo "**********Reference dictionary creation completed**********"
fi

# 3. BQSR
mkdir -p $out_dir/gatk/bqsr
echo "base quality score recalibrating..."
while IFS=',' read -r name input1 input2 condition sex; do
    bam=$out_dir/bam/$name.prcsd.bam
    $GATK BaseRecalibrator \
        -I $bam \
        -R $ref \
        --known-sites $dbsnp \
        --known-sites $mills \
        --known-sites $g1000 \
        -O $out_dir/gatk/bqsr/${name}_recal.table
    $GATK ApplyBQSR \
        -R $ref \
        -I $bam \
        --bqsr-recal-file $out_dir/gatk/bqsr/${name}_recal.table \
        -O $out_dir/gatk/bqsr/${name}_bqsr.bam
done < $sample_table

# 4-1. call somatic SNP/INDEL by Mutect2
if [ $dect_mode == "somatic" ]; then
    # get tumor and normal bam
    while IFS=',' read -r name input1 input2 condition sex; do
        if [ $condition == "tumor" ]; then 
            tumor_bam=$out_dir/gatk/bqsr/${name}_bqsr.bam
        else
            control_bam=$out_dir/gatk/bqsr/${name}_bqsr.bam
            control_name=$name
        fi
    done < $sample_table
    # run Mutect2
    mkdir -p $out_dir/gatk/vcf_interm
    mkdir -p $out_dir/gatk/else
    echo "**********Running GATK Mutect2...**********"
    $GATK Mutect2 \
        $interval_para \
        -R $ref \
        -I $tumor_bam \
        -I $control_bam \
        -normal $control_name \
        --germline-resource $germline_resource \
        -pon $pon_resource \
        -O $out_dir/gatk/vcf_interm/mutect2_raw.vcf.gz \
        --f1r2-tar-gz $out_dir/gatk/else/f1r2.tar.gz
    $GATK LearnReadOrientationModel \
        -I $out_dir/gatk/else/f1r2.tar.gz \
        -O $out_dir/gatk/else/read_orientation_model.tar.gz
    $GATK GetPileupSummaries \
        $interval_para \
        -I $tumor_bam \
        -V $common_biallelic \
        -L $common_biallelic \
        -O $out_dir/gatk/else/tumor_pileups.table
    $GATK GetPileupSummaries \
        $interval_para \
        -I $control_bam \
        -V $common_biallelic \
        -L $common_biallelic \
        -O $out_dir/gatk/else/control_pileups.table
    $GATK CalculateContamination \
        -I $out_dir/gatk/else/tumor_pileups.table \
        -matched $out_dir/gatk/else/control_pileups.table \
        -tumor-segmentation $out_dir/gatk/else/segments.table \
        -O $out_dir/gatk/else/contamination.table
    $GATK FilterMutectCalls \
        $interval_para \
        -R $ref \
        -V $out_dir/gatk/vcf_interm/mutect2_raw.vcf.gz \
        --contamination-table $out_dir/gatk/else/contamination.table \
        --tumor-segmentation $out_dir/gatk/else/segments.table \
        --ob-priors $out_dir/gatk/else/read_orientation_model.tar.gz \
        -O $out_dir/gatk/vcf_interm/mutect2.vcf.gz
    gatk_vcf=$out_dir/gatk/vcf_interm/mutect2.vcf.gz
    echo "**********GATK Mutect2 execution completed**********"
fi

# 3-2. call germline SNP/INDEL by HaplotypeCaller
if [ $dect_mode != "somatic" ]; then
    echo "**********Running GATK HaplotypeCaller...**********"
    mkdir -p $out_dir/gatk/else/
    mkdir -p $out_dir/gatk/gvcf/
    mkdir -p $out_dir/gatk/vcf_interm/
    while IFS=',' read -r name input1 input2 condition sex; do
        sample=$name
        bam=$out_dir/gatk/bqsr/${name}_bqsr.bam
        echo "**********$sample calling gvcf**********"
        $GATK HaplotypeCaller \
            $interval_para \
            -I $bam \
            -R $ref \
            -ERC GVCF \
            -O $out_dir/gatk/gvcf/${sample}.g.vcf.gz
        command_gvcf_input+="-V $out_dir/gatk/gvcf/${sample}.g.vcf.gz "
    done < $sample_table
    echo "**********combine gvcfs**********"
    $GATK CombineGVCFs \
        $interval_para \
        -R $ref \
        $command_gvcf_input \
        -O $out_dir/gatk/gvcf/combine_cohort.g.vcf.gz
    echo "**********genotype gvcfs**********"
    $GATK GenotypeGVCFs \
        $interval_para \
        -R $ref \
        -V $out_dir/gatk/gvcf/combine_cohort.g.vcf.gz \
        -O $out_dir/gatk/vcf_interm/HC_raw.vcf.gz
    # echo "**********Variant Quality Score Recalibration**********"
    sample_num=$( less $sample_table |wc -l )
    if [ $sample_num -ge 10 ]; then para_vqsr1="-an InbreedingCoff"; fi
    if [ $seq_method == "wgs" ]; then para_vqsr2="-an DP"; fi
    $GATK VariantRecalibrator \
        $interval_para \
        -R $ref \
        -V $out_dir/gatk/vcf_interm/HC_raw.vcf.gz \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $g1000 \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR $para_vqsr1 $para_vqsr2 \
        -mode SNP \
        -O $out_dir/gatk/else/snp.recal \
        --tranches-file $out_dir/gatk/else/snp.tranches
    $GATK ApplyVQSR \
        $interval_para \
        -R $ref \
        -V $out_dir/gatk/vcf_interm/HC_raw.vcf.gz \
        -O $out_dir/gatk/vcf_interm/vqsr_1.vcf.gz \
        --truth-sensitivity-filter-level 99.5 \
        -mode SNP \
        --recal-file $out_dir/gatk/else/snp.recal \
        --tranches-file $out_dir/gatk/else/snp.tranches
    $GATK VariantRecalibrator \
        $interval_para \
        -R $ref \
        -V $out_dir/gatk/vcf_interm/vqsr_1.vcf.gz \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
        -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum $para_vqsr $para_vqsr2 \
        -mode INDEL \
        -O $out_dir/gatk/else/indel.recal \
        --tranches-file $out_dir/gatk/else/indel.tranches \
        --max-gaussians 4
    $GATK ApplyVQSR \
        $interval_para \
        -R $ref \
        -V $out_dir/gatk/vcf_interm/vqsr_1.vcf.gz \
        -O $out_dir/gatk/vcf_interm/HC.vcf.gz \
        --truth-sensitivity-filter-level 99 \
        -mode INDEL \
        --recal-file $out_dir/gatk/else/indel.recal \
        --tranches-file $out_dir/gatk/else/indel.tranches
    gatk_vcf=$out_dir/gatk/vcf_interm/HC.vcf.gz
    echo "**********GATK HaplotypeCaller execution completed**********"
fi

# 4. extract
$BCFTOOLS view \
    -f PASS \
    -v snps \
    -o $out_dir/gatk/snv_gatk.vcf \
    $gatk_vcf
$BCFTOOLS view \
    -f PASS \
    -v indels \
    -o $out_dir/gatk/indel_gatk.vcf \
    $gatk_vcf
