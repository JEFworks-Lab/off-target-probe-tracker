# `otpc`: off-target-probe-checker

`otpc` is a simple python program that aligns probe sequences to transcript sequences to detect potential off-target probe activities.

## Installation

`otpc` is optimized for linux systems. We recommend that the users install the program in a conda environment as follows:

```
$ conda create --name otpc pip
$ conda activate otpc
$ conda install bioconda::mummer4
$ conda install samtools
$ conda install bioconda::gffread bioconda::bowtie2
$ pip install .
```

## Usage

`otpc` consists of three modules: `flip`, `track`, and `stat`. `flip` corrects probes that are aligning to the opposite strand of their intended target genes, by reverse complementing them. We assume probe sequences are designed in the same strand as their targets. The module requires the annotation for the target transcripts as well as their sequences. We recommend that the users use [gffread](https://github.com/gpertea/gffread) to extract processed transcript sequences from annotation GFF/GTF files (e.g., `$ gffread -w transcripts.fa -g genome.fa transcripts.gff`).

```
$ otpc -o out_dir flip -i probes.fa -a transcripts.gff -f transcripts.fa
```

This module outputs forward oriented probe sequences in a file called `fwd_oriented.fa`. 

`track` is the main module that aligns query probe sequences to any target transcriptome. We recommend that the users be conservative when choosing their transcriptome to avoid false positive prediction of off-target binding. `otpc` predict off-target binding by aligning query probes to target transcripts. By default, binding is predicted for only perfect matches (i.e., no indels, clips, or mismatches). See options for flags that allow for more lenient predictions that allow misalignments.

```
$ otpc -o out_dir track -q query.fa -t target.fa -a target.gff # query.fa most likely will be fwd_oriented.fa
```

This module outputs a CSV file containing the gene and transcript information to which each probe aligns. Each probe is also annotated with the number of genes it aligns to as well as the CIGAR strings for its alignments.

```
$ otpc -o out_dir stat -i probe2targets.tsv -q query.fa
```

This module summarizes `otpc` binding predictions. For each targeted gene, stat.summary.tsv file shows the number of probes and the genes those probes aligns to. Each pair of (target_gene, binding_gene), the module annotates number of alignments to the binding_gene and the corresponding number of probes (n of probes << n of alignmennts). 

## Arguments

See below for the complete list of arguments:
```
Usage: otpc [common_args] [MODULE] [args]

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
