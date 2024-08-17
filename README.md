

# hDNApipe
Streamlining human genome analysis and interpretation with an intuitive and user-friendly pipeline




## Table of Contents


## Introduction




## Installation
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

## Command-Line usage
4 module

## Input (sample information table)

## Output

## config

## Graphical user interface

## Example workflows

# Troubleshooting
If you encounter errors from hDNApipe, please report them on the issues page. Any bug reports, suggestions and general feedback would be highly welcome.


