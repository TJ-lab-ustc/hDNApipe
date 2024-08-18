# hDNApipe
Streamlining human genome analysis and interpretation with an intuitive and user-friendly pipeline

<img src="https://github.com/user-attachments/assets/1816e1c0-dcbb-4df2-b660-bfdb43ebe002" width="700px">

## Table of Contents
1. [Introduction](#Introduction)
2. [Setup](#Serup)
3. [Usage](#Usage)
4. [Config](#Config)
5. [Input](#Input)
6. [Output](#Output)
7. [Example Workflow](#Example)


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

### Command-Line
There are four components in hDNApipe:
  1. init - This is used for downloading necessary resource files when installing hDNApipe for the first time. However, they have already been downloaded in the docker image. There is no necessity to utilize this, but we retained this function just in case.
  2. ref - This is used for setting the reference genome. 
  3. var - The main module in hDNApipe. It is used to run the genomic analysis pipeline, including alignment, preprocessing, variant calling, and annotation.
  4. plot - Used to plot additional statistics and analysis graphs.
```
usage:  dnapipe <command> [options]

command: init   Download the files necessary for variants calling.
         ref    Select a reference genome from an existing file or automatically download it. 
                An automatic indexing will be performed if related files are not found.
         var    Run several variant calling pipelines from a clean fastq file or bam file. 
         plot   Plot additional statistics and analysis graphs.

options: 
        -h       show this message
```

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
This is utilized for setting the reference genome. Two modes are provided. One is downloading the recommended reference online and indexing it; another one is choosing from an existing fasta file and indexing it if no related files are found in its directory.

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
There are basic arguments must be declared in running. 
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



### Graphical user interface

## Config

## Input (sample information table)

## Output




## Example workflow



# Troubleshooting
If you encounter errors from hDNApipe, please report them on the [issues](https://github.com/TJ-lab-ustc/hDNApipe/issues) page. Any bug reports, suggestions and general feedback would be highly welcome.


