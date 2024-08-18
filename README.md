# hDNApipe
Streamlining human genome analysis and interpretation with an intuitive and user-friendly pipeline

<img src="https://github.com/user-attachments/assets/1816e1c0-dcbb-4df2-b660-bfdb43ebe002" width="700px">

## Table of Contents
1. [Introduction](#Introduction)
2. [Setup](#Setup)
3. [Usage](#Usage)

   * [Command-Line](#Command-Line)

   * [Graphical user interface](#Graphical-user-interface)
   
4. [Input](#Input)
5. [Output](#Output)
6. [Troubleshooting](#Troubleshooting)


## Introduction




## Setup
Due to the complexity of configuring the environment necessary for hDNApipe, as it involves numerous tools and dependencies, we have made hDNApipe into a Docker container image based on the Ubuntu operating system.

The docker image can be downloaded from xxx.

First, use git to download the latest development tree. 
Then, the annotation files should be downloaded separately as they are not included in the Docker image due to their large size, which totals approximately 47GB. It is recommended to run the bash script used to download the annotation files in a tmux window or a similar manner, as it may take a considerable amount of time to complete. Moreover, it is suggested to download in the main computer and bind mount a volume instead of downloading in the docker container.

This download step can be skipped if annoation function is not needed.
```
git clone https://github.com/TJ-lab-ustc/hDNApipe
cd hDNApipe
chmod +x ./dnapipe
bash download_annoatation.sh <out_dir>
```

## Usage
  docker

### Command-Line
There are four components in hDNApipe.

1. [dnapipe init](#dnapipe-init) for the first time downloading.

2. [dnapipe ref](#dnapipe-ref) for setting reference genome.

3. [dnapipe var](#dnapipe-var) for genomic analysis pipeline.

4. [dnapipe plot](#dnapipe-plot) for plotting.

#### dnapipe init
This is used for downloading necessary resource files when installing hDNApipe for the first time. However, they have already been downloaded in the docker image. There is no necessity to utilize this, but we retained this function just in case.

Download files used in variant calling: 
```
dnapipe init --var [-o output_dir]
```
Download files used in annotations, which may take lots of time:
```
dnapipe init --annot [-o output_dir]
```

#### dnapipe ref
This is utilized for setting the reference genome. Two modes are provided. One is downloading the recommended reference online and indexing it; another one is choosing from an existing fasta file and indexing it if no related files are found in its directory. In fact, the recommended reference is also included in the docker image.

To download the recommended hg38 reference online:
```
dnapipe ref --download <save_dir>
```
To set an existing file as the reference:
```
dnapipe ref --set <ref.fa>
```

#### dnapipe var
The main module in hDNApipe. It is used to run the genomic analysis pipeline, including alignment, preprocessing, variant calling, and annotation.
```
usage:   dnapipe var <core_options> [extra_options]
```
There are basic arguments that must be declared during running var module. There are some options for the detection mode, sequencing method, sequencing input file type, and variant types to call. A region file is required when analyzing WES or target data. Sample information table should be needed as [Input](#Input). 
```
Core options:
        --mode/-M           Detection mode: germline or somatic
        --seq-method/-s     Sequence method: wgs or wes (For the target, choose wes).
        --file-type/-f      Input sequencing file type: fastq or bam
        --variant/-v        Variant type to call: short (refers to snv and indel), cnv or sv.
                                (Use , to link, such as short,cnv,sv)
        --sample-info/-i    Provide sample informtion table.
        --region/-r         Region to detect variants on. Bed file. Only necessary for wes/target data.
```
And there are more parameters for advanced users to customize: 
```
Computing options:
        -t threads  INT     Threads employed for multi-thread tasks (default: 1/2 of the total)
        -m memory   INT     Set the maximum memory in G (default: 1/2 of the total)

Align and preprocess:
        --mapq/-Q   INT     Minimum map quality score for output (default: 30)
        --remove-dup/-D     Remove duplicates rather than just marking them. (default: disabled)
        --soft-clip/-Y      Use soft clipping for supplementary alignments instead of hard clipping. (default: disabled)
        --force-rg          Force the addition or replacement of read group information. (default: disabled)

Variants calling:
        --short   STR       Select SNP/Indel calling tools: strelka, deepvariant, gatk.  (default: strelka)
                                Note that deepvariant does NOT support somatic.
        --sv      STR       Select SV calling tools: manta, lumpy, delly (default: manta)
        --cnv     STR       Select CNV calling tools: cnvkit, delly (default: cnvkit)
        --bin     INT       CNV calling bin size in bp (default: 10000)

Annotation: (default: options not enabled)
        --annot             Use VEP to annotate SNP/INDEL vcf with default information and AnnotSV for SV/CNV                                                 )
    More VEP annotation information can be added:
        --symbol            Add Gene Symbol information.
        --exist             Check the existence of each variant.
        --pathogenicity     Predict the pathogenicity of SNV and INDEL by SIFT, PolyPhen2, and CADD respectively.
        --frequency         Add Population allele frequency.
        --clinvar           Add Clinvar database annotation.
    Reports:
        --report            Generate reports for annotation results.
    Else:
        --annot-all         Use all the above options in annotation.
```

#### dnapipe plot
Used to plot additional statistics and analysis graphs from vcf files. Due to the differences between VCFs from different variant calling tools, plots including pathogenicity, GO, KEGG, and PPI analysis are limited to short variants annotated by VEP with pathogenicity annotations, and the CNV VCF used for circos is only supported by CNVkit output vcf.
```
usage:   dnapipe plot <input_options> <plot_types> [extra_options]

Input options:
        --snp           Provide SNP vcf.
        --indel         Provide INDEL vcf.
        --cnv           Provide CNV vcf.
        --sv            Provide SV vcf.
        -o/--outdir     Specify the output directory. (default: the current working directory)

Plot types:
        --category      Plot the category.
        --length        Plot the length distribution.
        --pathogenicity Plot the pathogenicity prediction.
        --circos        Plot the Circos plot.
        --go, --kegg, --ppi
                        Plot the analysis result of GO/KEGG/PPI analysis.

Additional Options:
        --bin           Set the bin size for the Circos plot, in bp. (default: 5000000)
        --cadd_score    Set the score threshold for CADD. (default: 10)
        -h/--help       Show this help message and exit.
```

### Graphical user interface


## Input
The resources required for hDNApipe are declared in the config file, and the template is placed in [dnapipe.config]([url](https://github.com/TJ-lab-ustc/hDNApipe/blob/main/dnapipe.config)). If no modification of the docker container is made, the only change for the user is to complete three paths in the config: `dnapipe_dir`, `dir_vep_annot` and `annotsv_dir`. For example:
```
dnapipe_dir="/opt/hDNApipe"
dir_vep_annot="/opt/download/vep_annot"
annotsv_dir="/opt/download/AnnotSV_annotations"
```

The sample information table is required for `dnapipe var`. It should contain information including: sample name, path1, path2, condition and sex. The sample name and at least one input path are necessary. The sample name refers to the name used for adding the BAM read group and the prefix of most output files. The path is prepared for the sequencing file location.The condition is needed to make it clear which one is `tumor` and which one is `control` when running somatic analysis; otherwise, it is not needed. Sex is optional. Use ',' as seperator. For example, a full information table looks like:
```
sample1,/path/sample1_1.fastq.gz,/path/sample1_2.fastq.gz,tumor,male
sample2,/path/sample2_1.fastq.gz,/path/sample2_2.fastq.gz,normal,male
```
If there is some unnecessary information missing, do it like: 
```
sample,/path/sample.bam,,,
```



## Output





# Troubleshooting
If you encounter errors from hDNApipe, please report them on the [issues](https://github.com/TJ-lab-ustc/hDNApipe/issues) page. Any bug reports, suggestions and general feedback would be highly welcome.
