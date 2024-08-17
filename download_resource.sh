# download all needed data files for hDNApipe

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

# annotsv

annotsv_dir=$download_dir/AnnotSV_annotations

mkdir $annotsv_dir

cd $download_dir

git clone https://github.com/lgmgeo/AnnotSV.git

cd AnnotSV

make PREFIX=. install

make PREFIX=. install-human-annotation

mv share/AnnotSV/Annotations_Exomiser ..

mv share/AnnotSV/Annotations_Human ..

cd ..

rm -r AnnotSV

# vep plugin

vep_annot_dir=$download_dir/vep_annot

mkdir -p $vep_annot_dir/Plugins

wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/112/CADD.pm \
    -O $vep_annot_dir/Plugins/CADD.pm

# vep clinvar + cadd

mkdir -p $vep_annot_dir/files

wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz \
    -O $vep_annot_dir/files/clinvar.vcf.gz

wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi \
    -O $vep_annot_dir/files/clinvar.vcf.gz.tbi

wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz \
    -O $vep_annot_dir/files/gnomad.genomes.r4.0.indel.tsv.gz

wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi \
    -O $vep_annot_dir/files/gnomad.genomes.r4.0.indel.tsv.gz.tbi

# edit config (download dir)

sed -i "s|^download_dir=\".*\"|download_dir=\"$download_dir\"|" $config

sed -i "s|^annotsv_dir=\".*\"|annotsv_dir=\"$annotsv_dir\"|" $config

plugin_dir=$download_dir/vep_file

sed -i "s|^dir_plugins=\".*\"|dir_plugins=\"$plugin_dir\"|" $config

# vep annotation (take long time)

mkdir -p $vep_annot_dir/Cache

# wget -c --tries=0 --read-timeout=10 http://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz \
#     -O $vep_annot_dir/Cache/homo_sapiens_vep_112_GRCh38.tar.gz

# if [ $? -eq 0 ]; then
#     tar -zxvf $vep_annot_dir/Cache/homo_sapiens_vep_112_GRCh38.tar.gz
#     rm $vep_annot_dir/Cache/homo_sapiens_vep_112_GRCh38.tar.gz
# fi

