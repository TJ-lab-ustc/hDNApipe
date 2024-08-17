#!/bin/bash

download_dir=$1

dnapipe_dir=$2

config=$dnapipe_dir/dnapipe.config

db_dir=$download_dir/db

mkdir -p $db_dir

# delly

wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz \
    -O $db_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz

wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.fai \
    -O $db_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.fai

wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.gzi \
    -O $db_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.gzi

wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv \
    -O $db_dir/human.hg38.excl.tsv

# cnvkit

wget https://raw.githubusercontent.com/etal/cnvkit/master/data/refFlat_hg38.txt \
    -O $db_dir/refFlat_hg38.txt

# gatk

eval "$(conda shell.bash hook)"

conda activate $dnapipe_env_py3

pip install gsutil

bundle_dir=gs://gcp-public-data--broad-references/hg38/v0

somatic_dir=gs://gatk-best-practices/somatic-hg38

gsutil cp $bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz* $db_dir

gsutil cp $bundle_dir/1000G_omni2.5.hg38.vcf.gz* $db_dir 

gsutil cp $bundle_dir/hapmap_3.3.hg38.vcf.gz* $db_dir

gsutil cp $bundle_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz* $db_dir 

gsutil cp $bundle_dir/Homo_sapiens_assembly38.dbsnp138.vcf.gz* $db_dir

gsutil cp $somatic_dir/small_exac_common_3.hg38.vcf.gz* $db_dir

gsutil cp $somatic_dir/1000g_pon.hg38.vcf.gz* $db_dir

gsutil cp $somatic_dir/af-only-gnomad.hg38.vcf.gz* $db_dir


# edit config (download dir)

sed -i "s|^download_dir=\".*\"|download_dir=\"$download_dir\"|" $config
