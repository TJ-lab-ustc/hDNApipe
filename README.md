# hDNApipe
Streamlining human genome analysis and interpretation with an intuitive and user-friendly pipeline.

<img src="https://github.com/user-attachments/assets/1816e1c0-dcbb-4df2-b660-bfdb43ebe002" width="700px">

## Table of Contents
1. [Introduction](#Introduction)
2. [Setup](#Setup)
3. [Usage](#Usage)

   * [Command-Line](#Command-Line)

   * [Graphical user interface](#Graphical-user-interface)
   
4. [Input](#Input)
5. [Output](#Output)
6. [Usage example and details](#Usage-example-and-details)
7. [Troubleshooting](#Troubleshooting)


## Introduction
hDNApipe is a highly flexible end-to-end pipeline designed for the analysis and interpretation of human genomic sequencing data. This tool is capable of detecting a wide range of variant types in both germline and somatic contexts, including single nucleotide variants (SNVs), small insertions and deletions (INDELs), large structural variants (SVs), and specifically copy number variations (CNVs). It has a dual-mode operation through both command-line and graphical user interface (GUI), ensuring an accessible user experience.

## Setup
Given the intricate nature of setting up the environment required for hDNApipe, which entails a multitude of tools and dependencies, we have encapsulated hDNApipe into a Docker container image predicated on the Ubuntu operating system. The corresponding Dockerfile will be made available soon.

First, utilize the git command to retrieve the most recent development tree and build Docker image from Dockerfile. 
```
git clone https://github.com/TJ-lab-ustc/hDNApipe
cd hDNApipe
chmod +x ./dnapipe
docker build -t hDNApipe:v1.0 -f Dockerfile
```

Subsequently, obtain the files associated with the annotation function. This download step can be skipped if annotation function is not needed.

It should be noted that the annotation files, amounting to approximately 47GB in total, are not incorporated within the Docker image due to their substantial size. Hence, they need to be downloaded independently. It is advisable to execute the bash script designed for downloading these annotation files within a **tmux** window or an analogous environment, considering that the download process might consume a significant amount of time. Furthermore, it is recommended to conduct the download on the host machine and then mount a volume, as opposed to downloading directly within the docker container. 
```
bash download_annoatation.sh <out_dir>
```
Specify the location where the annotation information will be downloaded in ```<out_dir>```.

## Usage
Use `docker run -v` to to map into the container the folder of our tool `hDNApipe`, annotation files folder, and your sequencing data folder.  
```
docker run \
  -v /path/hDNApipe:/hDNApipe \
  -v /path/annotation_dir:/annotation_dir \
  -v /path/data:/input \
  -it hDNApipe /bin/bash
```
When running, change ```/path/hDNApipe/``` to the folder path of the hDNApipe downloaded via Git. Change ```/path/annotation_dir/``` to the folder path specified when downloading the annotation files. And ```/path/data``` should be set as the folder path which contains all the sequencing files.
*The `-v` option in the command is used for volume mounting. It allows you to map directories from your local machine to directories inside the Docker container. The part before the colon (:) refers to a directory path on your local machine, and the part after the colon indicates the corresponding directory path inside the Docker container. 

When aiming to enable access to the remote GUI, certain crucial additional settings need to be configured:
```
--net=host -e DISPLAY=:10.0 -v /$HOME/.Xauthority:/root/.Xauthority
```
The ```--net=host``` option allows the container to share the host's network namespace, which is essential for proper communication and access to the GUI. The ```-e DISPLAY=:10.0``` environment variable setting specifies the display server to which the GUI applications will connect. This ensures that the graphical output is directed to the correct display. The ```-v /$HOME/.Xauthority:/root/.Xauthority``` volume mount is used to share the X authority file between the host and the container, which is necessary to authenticate and authorize the container to access the host's X server, enabling a seamless and secure GUI experience when working with remote systems or containers.

After entering Docker, operations can be divided into two forms: command line and graphical user interface.

### Command-Line
There are three main subcommands in hDNApipe.

1. [dnapipe ref](#dnapipe-ref) for setting reference genome.

2. [dnapipe var](#dnapipe-var) for genomic analysis pipeline.

3. [dnapipe plot](#dnapipe-plot) for plotting.

#### dnapipe ref
This functionality is designed to configure the reference genome, with two available modes. The first mode entails downloading the recommended reference from the online repository and subsequently generating an index for it. The second mode involves selecting from an existing fasta file; if no relevant index files are detected within its directory, an indexing process will be initiated. 
We recommend using the recommended reference genome to circumvent the issue of incompatibility with certain individual tools. Notably, the recommended reference genome is already part of the docker image, so actually there's no need to download it under normal circumstances.

To set an existing file as the reference:
```
dnapipe ref --set <ref.fa>
```
```<ref.fa>```is used to specify the reference genome FASTA file.

To download the recommended hg38 reference online:
```
dnapipe ref --download <save_dir>
```
```<save_dir>```is used to specify the path when downloading the reference genome


#### dnapipe var
This is the main module in hDNApipe. It is used to run the genomic analysis pipeline, including alignment, preprocessing, variant calling, and annotation.
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
For example, run the following code to perform germline detection on WES fastq data, and the mutation types to be detected are specified as short variants (SNV and INDEL) and SV. All the results will be placed into the ```/output_dir```.
```
./dnapipe var \
  	--mode germline --file-type fastq --seq-method wes \
  	--region /input/idt_capture_novogene_no_alt.bed \
  	-i /input/sample.txt \
  	--variant short,sv \
  	-o /output_dir
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
        --annot             Use VEP to annotate SNP/INDEL vcf with default information and AnnotSV for SV/CNV                                           
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
For example, modify the number of threads by adding the ```-t``` parameter, and ```--short``` can be used to specify the caller for detecting short variants:
```
./dnapipe var \
  	--mode germline --file-type fastq --seq-method wes \
  	--region /input/idt_capture_novogene_no_alt.bed \
  	-i /input/sample.txt \
  	--variant short,sv \
  	-o /output_dir \
  	-t 20 --short deepvariant
```

#### dnapipe plot
Used to plot additional statistics and analysis graphs from vcf files. Due to the differences between VCFs from different variant calling tools, plots including pathogenicity, GO, KEGG, and PPI analysis are limited to short variants annotated by VEP with pathogenicity annotations, and the CNV VCF used for circos is only supported by CNVkit output vcf. There is an expectation for the expansion of its compatibility down the road.
```
usage:   dnapipe plot <input_options> <plot_types> [extra_options]

Input options:
        --snp           Provide SNV vcf.
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
For example, the following can be used to plot the classification information of SNV.

```dnapipe plot --snp /path/snv.vcf --category -o /path/outdir```

Here are some example plots for the visualization function:
![image](https://github.com/user-attachments/assets/32085cf6-2edb-45c7-8f3e-edcd05de0084)


### Graphical user interface
Pop up the GUI window by entering the command 'python dnapipe.py', and subsequent operations do not require the command line. Its main interface is shown in the figure below：
![image](https://github.com/user-attachments/assets/4f0769ca-bac1-4ca5-9f22-db1425c04738)

Upon the initial usage, the reference genome can be configured via the "Reference" button located within the "Initialization" tab. Then, parameters can be set on the two options pages. The basic options page contains the essential parameters, while the advanced options page has the optional ones. The options in the GUI correspond to the arguments in the command-line. After completion, just click the RUN button to start the operation.

Details on GUI usage and examples can be found in the [manual instruction](https://github.com/TJ-lab-ustc/hDNApipe/blob/main/manual.pdf).

## Input
The resources required for hDNApipe are declared in the config file, and the template is placed in [dnapipe.config]([url](https://github.com/TJ-lab-ustc/hDNApipe/blob/main/dnapipe.config)). If no modification of the docker container is made, the only change for the user is to complete two paths in the config: `dnapipe_dir`, `dir_vep_annot` and `annotsv_dir`. For example:
```
dnapipe_dir="/hDNApipe"                # hDNApipe home directory
annotation_dir="/annotation_dir"       # annotation files download directory
```
Note that both of these require the location within the Docker container, that is, the location after the colon when running the ```docker -v``` command. 

The sample information table is required for `dnapipe var`. It should contain information including: sample name, path1, path2, condition and sex. The sample name and at least one input path are necessary. The sample name refers to the name used for adding the BAM read group and the prefix of most output files. The path is prepared for the sequencing file location.The condition is needed to make it clear which one is `tumor` and which one is `control` when running somatic analysis; otherwise, it is not needed. Sex is optional. Use ',' as seperator. For example, a full information table looks like:
```
sample1,/path/sample1_1.fastq.gz,/path/sample1_2.fastq.gz,tumor,male
sample2,/path/sample2_1.fastq.gz,/path/sample2_2.fastq.gz,control,male
```
If there is some unnecessary information missing, do it like: 
```
sample,/path/sample.bam,,,
```


## Output
The ultimate outputs consist of annotated variant call format files, visualization figures and reports. However, hDNApipe also saves intermediate files, which facilitates users to review and check. The structure of the output files is as follows:
![image](https://github.com/user-attachments/assets/b376ee85-88bd-43a3-ad45-18f5fdb776c1)

## Usage example and details 
For practical operation examples, please refer to the [manual instruction](https://github.com/TJ-lab-ustc/hDNApipe/blob/main/manual.pdf).

# Troubleshooting
If you encounter errors from hDNApipe, please report them on the [issues](https://github.com/TJ-lab-ustc/hDNApipe/issues) page. Any bug reports, suggestions and general feedback would be highly welcome.
