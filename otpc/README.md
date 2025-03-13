# `otpc`: off-target-probe-checker

`otpc` is a simple python program that aligns probe sequences to transcript sequences to detect potential off-target probe activities.

## Installation

`otpc` is optimized for linux systems. We recommend that the users install the program in a conda environment.

```
# substitute otpc with another env name
$ conda create --name otpc pip -y
$ conda activate otpc
$ conda install samtools bioconda::bowtie2 bioconda::gffread -y
$ pip install .
```
copy-and-paste above commands sequentially to set up the environment.

## Usage

Below is a sample command:

```
# --fwd flag enables automatic reverse complementation of reversely oriented probes
otpc -q probes.fa -t transcripts.fa -a transcripts.gff --fwd -o out [--nucmer]
```

TODO:
- install nucmer on pip install... this took me a bit to figure out lol
