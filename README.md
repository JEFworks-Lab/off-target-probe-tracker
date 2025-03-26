# `opt`: <ins>O</ins>ff-target <ins>P</ins>robe <ins>T</ins>racker

`opt` is a simple python program that aligns probe sequences to transcript sequences to detect potential off-target probe activities.

## Installation

`opt` is optimized for linux systems. We recommend that the users install the program in a conda environment as follows:

```
conda create --name opt pip python=3.9
conda activate opt
conda install bioconda::mummer4
conda install bioconda::gffread
conda install bioconda::samtools
conda install -c bioconda -c conda-forge bowtie2
git clone git@github.com:JEFworks-Lab/off-target-probe-tracker.git
cd off-target-probe-tracker/
pip install .
```

```
Wget https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.g
tar -xvzf mummer-4.0.1.tar.gz
Cd mummer-4.0.1
./configure --prefix=$(pwd)
Make
Make install
export PATH=$PATH:$(pwd)
```

## Usage

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

## Note

The target gene name and ID within the query.fa is expected to be in the following format:

`>gene_id|gene_name|accession`

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
          TODO
      --keep-dot
          TODO
      --force
          TODO
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
