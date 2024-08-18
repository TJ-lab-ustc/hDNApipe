#!/bin/bash
set -e

# function: call germline or somatic SV by lumpy
# usuage: script_dir/lumpy.sh $dect_mode $sample_table $out_dir $threads
# note: only work on WGS data; exclude specific regions such as telomere; no region parameter

# 1. Import variables 
dect_mode=$1
sample_table=$2
out_dir=$3
threads=$4
script_dir=$(readlink -f `dirname $0`)
dnapipe_dir=$(dirname $script_dir)
config="$dnapipe_dir/dnapipe.config"
source $config

#  2. Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py2

# 3. Pre-process
while IFS=',' read -r name input1 input2 condition sex; do
    echo "Preprocessing for LUMPY ..."
    bam=$out_dir/bam/$name.prcsd.bam
    samtools view -b -F 1294 -@ $threads $bam -o $out_dir/bam/$name.disc.bam
    samtools view -h -@ $threads $bam | \
        extractSplitReads_BwaMem -i stdin | \
        samtools view -Sb -@ $threads - \
        > $out_dir/bam/$name.split.bam
done < $sample_table

# 4. Run LUMPY
mkdir -p $out_dir/lumpy/vcf_interm
# 4-1 somatic
if [ $dect_mode == "somatic" ]; then
    bams=""
    # get tumor and normal bam
    while IFS=',' read -r name input1 input2 condition sex; do
        if [ $condition == "tumor" ]; then 
            tumor_bam=$out_dir/bam/$name.prcsd.bam
            tumor_split=$out_dir/bam/$name.split.bam
            tumor_disc=$out_dir/bam/$name.disc.bam
        else
            control_bam=$out_dir/bam/$name.prcsd.bam
            control_split=$out_dir/bam/$name.split.bam
            control_disc=$out_dir/bam/$name.disc.bam
        fi
        bam=$out_dir/bam/$name.prcsd.bam
        bams+="$bam,"
        split=$out_dir/bam/$name.split.bam
        splits+="$split,"
    done < $sample_table
    # run
    lumpyexpress \
        -B $tumor_bam,$control_bam \
        -S $tumor_split,$control_split \
        -D $tumor_disc,$control_disc \
        -o $out_dir/lumpy/vcf_interm/tumor_normal.vcf
    lumpy_vcf=$out_dir/lumpy/vcf_interm/tumor_normal.vcf
    bams="$tumor_bam,$control_bam"
# 4-2 germline
else
    bams=""
    splits=""
    discs=""
    while IFS=',' read -r name input1 input2 condition sex; do
        bam=$out_dir/bam/$name.prcsd.bam
        split=$out_dir/bam/$name.split.bam
        disc=$out_dir/bam/$name.disc.bam
        bams+="$bam,"
        splits+="$split,"
        discs+="$disc,"
    done < $sample_table
    lumpyexpress \
        -B $bams \
        -S $splits \
        -D $discs \
        -o $out_dir/lumpy/vcf_interm/samples.vcf
    lumpy_vcf=$out_dir/lumpy/vcf_interm/samples.vcf
    bams=${bams%,}
fi

# 5. Genotype
echo "Genotyping for LUMPY output ..."
awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' $ref.fai > $out_dir/lumpy/contigs.txt
cat $out_dir/lumpy/contigs.txt | sed '3r /dev/stdin' $lumpy_vcf > $out_dir/lumpy/samples_header.vcf
lumpy_vcf=$out_dir/lumpy/samples_header.vcf
svtyper \
    -B $bams \
    -i $lumpy_vcf \
    -o $out_dir/lumpy/vcf_interm/genotyped.vcf
rm $out_dir/lumpy/contigs.txt
rm $out_dir/lumpy/samples_header.vcf

# 6. Change conda environment
eval "$(conda shell.bash hook)"
conda activate $dnapipe_env_py3

# 7. Manual filter
bgzip $out_dir/lumpy/vcf_interm/genotyped.vcf
bcftools sort \
    $out_dir/lumpy/vcf_interm/genotyped.vcf.gz \
    -O z \
    -o $out_dir/lumpy/vcf_interm/genotyped.sort.vcf.gz
if [ $dect_mode == "somatic" ]; then
    bcftools view \
        -i '(GT[1]="mis" || GT[1]="0/0") && (GT[0]!="mis" && GT[0]!="0/0")' \
        $out_dir/lumpy/vcf_interm/genotyped.sort.vcf.gz \
        -o $out_dir/lumpy/sv_lumpy.vcf
else
    bcftools view \
        -e '( GT[*]="mis" || GT[*]="0/0")' \
        $out_dir/lumpy/vcf_interm/genotyped.sort.vcf.gz \
        -o $out_dir/lumpy/sv_lumpy.vcf
fi
