# Conda
FROM continuumio/miniconda3 as conda_setup
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge
RUN conda create -n py27 \
                    python=2.7.15 \
                    bioconda::strelka \
                    bioconda::manta \
                    bioconda::lumpy-sv \
                    bioconda::svtyper \
    && conda clean -a

RUN conda create -n py38 \
                    python=3.8.10 \
                    bioconda::cnvkit \
                    bioconda::annotsv \
		    bioconda::ensembl-vep \
		    bwa samtools sambamba \
		    bedtools bcftools \
    && conda clean -a


# GATK
FROM broadinstitute/gatk:4.5.0.0 as gatk

# Use google deepvariant 1.6.1
FROM google/deepvariant:1.6.1
COPY --from=conda_setup /opt/conda /opt/conda
RUN conda init && \
    conda config --set auto_activate_base false
COPY --from=gatk /gatk /gatk

# reference and db
RUN mkdir -p /opt/download/reference
ENV hg38_dir=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
RUN wget $hg38_dir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai -O /opt/download/reference/GRCh38_no_alt_analysis_set.fna.fai && \
    wget $hg38_dir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz && \
    wget $hg38_dir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz -O /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz && \
    gunzip -c /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > /opt/download/reference/GRCh38_no_alt_analysis_set.fna && \
    tar -xzf /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz && \
    mv /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.amb $ref_dir/GRCh38_no_alt_analysis_set.fna.amb && \
    mv /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.ann $ref_dir/GRCh38_no_alt_analysis_set.fna.ann && \
    mv /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwt $ref_dir/GRCh38_no_alt_analysis_set.fna.bwt && \
    mv /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.pac $ref_dir/GRCh38_no_alt_analysis_set.fna.pac && \
    mv /opt/download/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sa $ref_dir/GRCh38_no_alt_analysis_set.fna.sa

# base
RUN apt update && \
    apt install -y \
    openjdk-17-jdk \
    vim apt-utils \
    bc python3-tk \
    xauth xorg \
    && apt clean \
    && rm -rf /var/lib/apt/lists/*

# R packages
RUN conda install -n py38 r-base=4.3.3 delly=1.2.6 \
		    bioconda::glnexus
RUN conda install -n py38 \
                    r-matrix r-mass r-ggraph r-ggplot2 \
                    r-reshape2 r-scales r-dplyr r-gtable \
                    r-tidyr r-gridBase r-BiocManager r-magrittr \
                    r-stringr r-circlize r-httr2 \
                    bioconda::bioconductor-complexheatmap \
                    bioconda::bioconductor-clusterprofiler \
                    bioconda::bioconductor-stringdb \
                    bioconda::bioconductor-dnacopy \
                    bioconda::bioconductor-org.hs.eg.db \
    && conda clean -a

# delly Rscript
RUN mkdir -p /opt/delly/R
COPY ./deps/cnv.R /opt/delly/R
COPY ./deps/rd.R /opt/delly/R

# knotAnnotSV 
COPY ./deps/knotAnnotSV.zip /tmp
RUN unzip /tmp/knotAnnotSV.zip -d /opt/ && \
    mv /opt/knotAnnotSV-master /opt/knotAnnotSV && \
    rm /tmp/knotAnnotSV.zip
RUN conda run -n py38 cpan YAML::XS && \
    conda run -n py38 cpan Sort::Key::Natural

# variantconvert
COPY ./deps/variantconvert.zip /tmp
RUN unzip /tmp/variantconvert.zip -d /opt/ && \
    mv /opt/variantconvert-master /opt/variantconvert && \
    rm /tmp/variantconvert.zip
RUN cd /opt/variantconvert && \
    conda run -n py38 pip install -e .
RUN conda run -n py38 variantconvert init

# deps(delly, cnvkit)
RUN mkdir -p /opt/download/db
RUN wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz -O /opt/download/db/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz && \
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.fai -O /opt/download/db/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.fai && \
    wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.gzi -O /opt/download/db/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz.gzi && \
    wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv -O /opt/download/db/human.hg38.excl.tsv && \
    wget https://raw.githubusercontent.com/etal/cnvkit/master/data/refFlat_hg38.txt -O /opt/download/db/refFlat_hg38.txt

# deps(gatk)
RUN source activate py38
ENV bundle_dir=gs://gcp-public-data--broad-references/hg38/v0
ENV somatic_dir=gs://gatk-best-practices/somatic-hg38
RUN pip install gsutil && \
    gsutil cp $bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz* /opt/download/db/ &&\
    gsutil cp $bundle_dir/1000G_omni2.5.hg38.vcf.gz* /opt/download/db/ &&\
    gsutil cp $bundle_dir/hapmap_3.3.hg38.vcf.gz* /opt/download/db/ &&\
    gsutil cp $bundle_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz* /opt/download/db/ &&\
    gsutil cp $bundle_dir/Homo_sapiens_assembly38.dbsnp138.vcf.gz* /opt/download/db/ &&\
    gsutil cp $somatic_dir/small_exac_common_3.hg38.vcf.gz* /opt/download/db/ &&\
    gsutil cp $somatic_dir/1000g_pon.hg38.vcf.gz* /opt/download/db/ &&\
    gsutil cp $somatic_dir/af-only-gnomad.hg38.vcf.gz* /opt/download/db/

ENV PATH /gatk:$PATH

WORKDIR /opt

CMD ["/bin/bash"]
