#!/bin/bash
set -e

# usage: bash download_annotation.sh <download_dir>

if [ -z "$1" ]; then

    download_dir=$PWD

else

    download_dir=$(readlink -f $1)

fi

dnapipe_dir=`basename $0`

dnapipe_dir=$(readlink -f $dnapipe_dir)

config=$dnapipe_dir/dnapipe.config

# annotsv

annotsv_dir=$download_dir/AnnotSV_annotations

mkdir -p $annotsv_dir

cd $download_dir

git clone https://github.com/lgmgeo/AnnotSV.git

cd AnnotSV

make PREFIX=. install

make PREFIX=. install-human-annotation

mv share/AnnotSV/Annotations_Exomiser ..

mv share/AnnotSV/Annotations_Human ..

cd ..

rm -r AnnotSV

# edit config

sed -i "s|^annotsv_dir=\".*\"|annotsv_dir=\"$annotsv_dir\"|" $config

# vep plugin

vep_dir=$download_dir/vep_annot

# mkdir -p $vep_dir/Plugins

wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/112/CADD.pm \
    -O $vep_dir/Plugins/CADD.pm

# vep clinvar + cadd

mkdir -p $vep_dir/files

wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz \
    -O $vep_dir/files/clinvar.vcf.gz

wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi \
    -O $vep_dir/files/clinvar.vcf.gz.tbi

wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz \
    -O $vep_dir/files/gnomad.genomes.r4.0.indel.tsv.gz

wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi \
    -O $vep_dir/files/gnomad.genomes.r4.0.indel.tsv.gz.tbi

# edit config

sed -i "s|^vep_dir=\".*\"|vep_dir=\"$vep_dir\"|" $config

# vep annotation (take long time)

mkdir -p $vep_dir/Cache

wget -c --tries=0 --read-timeout=10 http://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz \
    -O $vep_dir/Cache/homo_sapiens_vep_112_GRCh38.tar.gz

tar -zxvf $vep_dir/Cache/homo_sapiens_vep_112_GRCh38.tar.gz

rm $vep_dir/Cache/homo_sapiens_vep_112_GRCh38.tar.gz
