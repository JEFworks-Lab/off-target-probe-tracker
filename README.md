# `opt`: <ins>O</ins>ff-target <ins>P</ins>robe <ins>T</ins>racker

`opt` is a simple python program that aligns probe sequences to transcript sequences to detect potential off-target probe activity.


## Installation

`opt` has been tested on Linux and Mac systems.


### Linux Installation

You will need to install the following packages and this repo. We recommend that the users install them in a new conda environment as follows:

```
conda create --name opt pip python=3.9
conda activate opt
conda config --add channels bioconda
conda config --add channels conda-forge
conda install gffread bowtie2 samtools mummer4
git clone git@github.com:JEFworks/off-target-probe-tracker.git
cd off-target-probe-tracker/
pip install .
```
Please check mummer4 version == 4.0.1


### Mac Installation

You will need to install the following packages and this repo. We recommend that the users install them in a new conda environment as follows:

```
conda create --name opt pip python=3.9
conda activate opt
conda config --add channels bioconda
conda config --add channels conda-forge
conda install gffread bowtie2 samtools
git clone git@github.com:JEFworks/off-target-probe-tracker.git
cd off-target-probe-tracker/
pip install .
```

To install mummer4 on Mac, you will need to use Brew rather than conda. Note that this will install it on your machine and not within the conda environment. To install mummer4 on Mac, use the following commands:

```
brew install autoconf automake libtool md5sha1sum
gem install yaggo
brew install mummer
```
These instructions can be found on the mummer [installation.md](https://github.com/mummer4/mummer/blob/master/INSTALL.md). Please check mummer4 version == 4.0.1


## Arguments

See below for the complete list of arguments:
```
Usage: opt [common_args] [MODULE] [args]

*common_args
      -o, --out-dir
          output directory (REQUIRED)
      -p, --threads
          number of threads
      --bam
          store alignment files as BAM instead of SAM
      -b
          binary path for aligners (bowtie2 or mummer)
      --gtf
          input annotation is in GTF format not GFF
      -l, --min-exact-match
          minimum exact match for mummer alignments
      --schema
          When loading an annotation file, the following five keys must be specified to
            define the schema used. These keys help extract essential transcript and
            gene information from the GTF/GFF file:
            1. feature type (3rd col) used for transcript entries
            2. transcript ID attribute (contained in 9th col)
            3. parent attribute for transcripts (contained in 9th col)
            4. gene name attribute (contained in 9th col)
            5. transcript type attribute (contained in 9th col)

            NOTE: annotations vary greatly in formats, so if you need assistance with
            determining which schema is appropriate, please open a git issue.
      --keep-dot
          TODO
      --force
          prevents OPT from skipping steps based on the existence of specific files from
            previous run (useful when the prev run was unsuccessful)
      --skip-index
          skip bowtie2 index building step
      
*flip args / options:
      -i, --in-file
          probe sequences fasta (REQUIRED)
      -a, --src-annotation
          source transcripts annotation (REQUIRED)
      -f, --src-fasta
          source transcript sequences fasta (REQUIRED)

*track args / options:
      -q, --query
          query probe sequences fasta (REQUIRED)
      -t, --target
          target transcript sequences fasta (REQUIRED)
      -a, --annotation
          target transcript annotation (REQUIRED)
      -1, --one-mismatch
          allow upto 1 mismatch
      -pl, --pad-length
          length of the pad where mis-alignment is allowed

*stat args / options:
      -i, --in-file
          track module results file (i.e., probe2targets.csv) (REQUIRED)
      -q, --query
          query probe sequences fasta (REQUIRED)
      --exclude-pseudo
          exclude pseudogenes when counting off-target probes and affected genes
      --pc-only
          only include protein coding genes
      -s, --syn-file
          gene synonyms CSV file with 2 columns
```


## Usage

There is a full example located in the [example.ipynb](https://github.com/JEFworks-Lab/off-target-probe-tracker/blob/main/example.ipynb) file. Below briefly describes what each module does.

`opt` consists of three modules: `flip`, `track`, and `stat`. 

`flip` corrects probes that are aligning to the opposite strand of their intended target genes by reverse complementing them. We assume probe sequences are designed in the same strand as their targets. The module requires the annotation for the target transcripts as well as their sequences. We recommend that the users use [gffread](https://github.com/gpertea/gffread) to extract processed transcript sequences from annotation GFF/GTF files (e.g., `$ gffread -w transcripts.fa -g genome.fa transcripts.gff`).

```
opt -o out_dir flip -i probes.fa -a transcripts.gff -f transcripts.fa
```

This module outputs forward oriented probe sequences in a file called `fwd_oriented.fa`. 

`track` is the main module that aligns query probe sequences to any target transcriptome. We recommend that the users be mindful of which target transcriptome they are using during this prediction step. `opt` predicts off-target binding by aligning query probes to target transcripts. By default, binding is predicted for only perfect matches (i.e., no indels, clips, or mismatches). See options for flags that allow for more lenient predictions that allow for misalignments.

Note that query.fa most likely will be fwd_oriented.fa

```
opt -o out_dir track -q query.fa -t target.fa -a target.gff
```

This module outputs a CSV file containing the gene and transcript information to which each probe aligns in a file called `probe2targets.tsv`. Each probe is also annotated with the number of genes it aligns to as well as the CIGAR strings for its alignments.

`stat` will summarize `opt` binding predictions.

```
opt -o out_dir stat -i probe2targets.tsv -q query.fa
```

For each targeted gene, the `stat.summary.tsv` file shows the number of probes and the genes those probes aligns to. For each pair of (target_gene, binding_gene), the module annotates number of alignments to the binding_gene and the corresponding number of probes (n of probes << n of alignmennts). Finally, the `collapsed_summary.tsv` file shows the target gene, number of probes, genes that the probes aligned to, number of alignments, and number of probes aligned to each gene in column 3 (similar to what is shown in Table 1 of our paper).


## Notes

### Probe ID format

The target gene name and ID within the query.fa is expected to be in the following format:

`>gene_id|gene_name|accession`


### mummer4 installation

It's important the mummer4 version is >= 4.0.1. If not, you can compile and install the latest release of mummer4 available [here](https://github.com/mummer4/mummer/releases). To compile and install mummer4:

if you've downloaded a newer release, replace 4.0.1 with the correct version number
```
wget https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz
tar -xvzf mummer-4.0.1.tar.gz
cd mummer-4.0.1
./configure --prefix=$(pwd)
make
make install
export PATH=$PATH:$(pwd)
```

To check if you've successfully installed MUMmer4, try running:

```
mummer -h
```

You should see the mummer help manual outputted in the terminal.

Note that every time you open a new kernel or shell session, you'll need to repeat the `EXPORT` command. To avoid it, you can add `export PATH=$PATH:$(pwd)` to your kernel / shell config file (e.g., `~/.bashrc`).

Similarly, if samtools is not installing through conda, we recommend that you compile and install it.
