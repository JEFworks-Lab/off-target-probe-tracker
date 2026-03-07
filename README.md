# OPT — Off-target Probe Tracker

OPT identifies off-target binding of spatial transcriptomics probes (Xenium, MERSCOPE, CosMx, and others) against a reference transcriptome using nucleotide alignment (nucmer or Bowtie2). It flags probes that align to genes other than their intended target and summarizes off-target activity by gene and transcript biotype.

**Citation:** Hallinan et al., *eLife* 2025. https://elifesciences.org/reviewed-preprints/107070


## Quick Start

OPT can be used via a **point-and-click web interface** (recommended) or directly from the **command line**.

---

## Installation

OPT has been tested on Linux and macOS. The fastest way to install is with the provided script:

```bash
git clone https://github.com/JEFworks-Lab/off-target-probe-tracker.git
cd off-target-probe-tracker
bash install.sh
```

This script will:
1. Create a conda environment named `opt` from `environment.yml`
2. Install mummer4 (via conda on Linux, via Homebrew on macOS)
3. Install the `opt` Python package

### Manual Installation

#### Linux

```bash
conda create --name opt pip python=3.9
conda activate opt
conda config --add channels bioconda
conda config --add channels conda-forge
conda install gffread bowtie2 samtools mummer4
git clone https://github.com/JEFworks-Lab/off-target-probe-tracker.git
cd off-target-probe-tracker
pip install .
```

#### macOS

```bash
conda create --name opt pip python=3.9
conda activate opt
conda config --add channels bioconda
conda config --add channels conda-forge
conda install gffread bowtie2 samtools
git clone https://github.com/JEFworks-Lab/off-target-probe-tracker.git
cd off-target-probe-tracker
pip install .
```

mummer4 must be installed via Homebrew on macOS (conda does not support it):

```bash
brew install autoconf automake libtool md5sha1sum
gem install yaggo
brew install mummer
```

> **Note:** mummer version >= 4.0.1 is required. Run `mummer -h` to confirm a successful install.

---

## Web Interface (Streamlit App)

The easiest way to use OPT is through the interactive web app:

```bash
conda activate opt
streamlit run app.py
```

Then open `http://localhost:8501` in your browser.


### App walkthrough

1. **Run Configuration** — set the output directory and number of threads.
2. **Input Files** — provide paths to your probe FASTA, target transcript FASTA, and annotation GFF/GTF. Use the Browse buttons or type paths directly. An optional gene synonyms CSV can be provided to remap gene names.
3. **Analysis Options:**
   - **Pad length** — number of bases at each probe end where mismatches are tolerated (0 = strict perfect match required in the core region).
   - **Max mismatches anywhere** — allow up to N mismatches anywhere across the full probe sequence. Can be combined with pad length: when both are set, both conditions must be satisfied.
   - **Exclude pseudogenes / Protein-coding only** — filter the off-target summary by biotype.
4. Click **Run OPT** to run all three modules (flip → track → stat) and view results in the dashboard below.

The results dashboard shows:
- Metric cards: genes with off-target binding, probes with off-target binding, protein-coding off-targets
- Gene-level off-target table (one row per target gene → off-target gene pair) with biotype badges and CIGAR strings, sortable and filterable by biotype
- Probe-level detail table (expandable)
- Download buttons for all key output files

---

## Command-Line Interface

OPT consists of three modules — `flip`, `track`, `stat` — plus an `all` module that runs all three in sequence.

### Run all modules at once (recommended)

```bash
opt -o out_dir all -q probes.fa -t transcripts.fa -a transcripts.gff
```

### Common arguments (apply to all modules)

| Argument | Description |
|---|---|
| `-o`, `--out-dir` | Output directory **(required)** |
| `-p`, `--threads` | Number of threads (default: 1) |
| `--bam` | Store alignments as BAM instead of SAM |
| `-l`, `--min-exact-match` | Minimum exact match length for nucmer (default: 20) |
| `--schema` | Comma-separated list of 5 GFF/GTF schema fields (see below) |
| `--keep-dot` | Keep version suffixes in gene IDs (e.g. ENSG00000.1) |
| `--force` | Recompute all steps, ignoring any cached results |
| `--skip-index` | Skip Bowtie2 index build step |

### `flip` — correct probe strand orientation

```bash
opt -o out_dir flip -q probes.fa -t transcripts.fa -a transcripts.gff
```

Probes are expected to be on the same strand as their target gene. `flip` detects probes that align to the reverse complement of their target and flips them. Output: `fwd_oriented.fa`.

### `track` — align probes and detect off-target binding

```bash
opt -o out_dir track -q fwd_oriented.fa -t transcripts.fa -a transcripts.gff
```

| Argument | Description |
|---|---|
| `-q`, `--query` | Query probe FASTA **(required)** |
| `-t`, `--target` | Target transcript FASTA **(required)** |
| `-a`, `--annotation` | Annotation GFF/GTF **(required)** |
| `-pl`, `--pad-length` | Tolerate mismatches in the terminal N bases of each probe end |
| `-mm`, `--max-mismatches` | Allow up to N mismatches anywhere in the full probe (-1 = disabled) |
| `-1`, `--one-mismatch` | Allow up to 1 mismatch using mummer exact-match extension |

Output: `probe2targets.tsv` (all probes) and `probe2targets_offtargets.tsv` (probes mapping to >1 gene).

### `stat` — summarize off-target predictions

```bash
opt -o out_dir stat -i probe2targets.tsv -q probes.fa
```

| Argument | Description |
|---|---|
| `-i`, `--in-file` | `probe2targets.tsv` from the track module **(required)** |
| `-q`, `--query` | Query probe FASTA **(required)** |
| `--exclude-pseudo` | Exclude pseudogenes from off-target counts |
| `--pc-only` | Count only protein-coding genes as off-targets |
| `-s`, `--syn-file` | Gene synonyms CSV (two columns: old name, new name) |

---

## Input File Formats

### Probe FASTA

Headers must follow this format:

```
>gene_id|gene_name|accession
```

Example:
```
>ENSG00000170458|CD14|22f9405
ATCGATCGATCGATCGATCG...
```

### Target transcript FASTA

Standard nucleotide FASTA of transcript sequences. We recommend extracting these with [gffread](https://github.com/gpertea/gffread):

```bash
gffread -w transcripts.fa -g genome.fa annotation.gff
```

### Annotation GFF/GTF

Standard GFF3 or GTF format, optionally gzip-compressed (`.gff.gz`, `.gtf.gz`). GENCODE, RefSeq, and CHESS formats are all supported. For non-standard annotation formats, use `--schema` to specify the correct field names (see below).

### Gene Synonyms CSV (optional)

Two-column CSV mapping old gene names to new names. No header required:

```
WARS,WARS1
CARS,CARS1
```

---

## GFF/GTF Schema

The `--schema` argument specifies five comma-separated field names used to parse the annotation:

```
transcript,ID,Parent,gene_name,transcript_type
```

| Position | Description | GENCODE GFF | RefSeq GFF | GTF |
|---|---|---|---|---|
| 1 | Feature type (col 3) | `transcript` | `transcript` | `transcript` |
| 2 | Transcript ID attribute | `ID` | `ID` | `transcript_id` |
| 3 | Parent (gene) attribute | `Parent` | `Parent` | `gene_id` |
| 4 | Gene name attribute | `gene_name` | `gene` | `gene_name` |
| 5 | Transcript type attribute | `transcript_type` | `transcript_biotype` | `transcript_type` |

If you are unsure which schema to use, open a GitHub issue.

---

## Output Files

| File | Description |
|---|---|
| `fwd_oriented.fa` | Strand-corrected probe sequences (from flip) |
| `flip_t2g.csv` | Transcript-to-gene map built during flip |
| `probe2targets.tsv` | All probe alignments with gene and CIGAR info |
| `probe2targets_offtargets.tsv` | Probes mapping to more than one gene |
| `collapsed_summary.tsv` | Per-gene summary of all probe alignments |
| `collapsed_summary_offtargets.tsv` | Per-gene summary for off-target genes only |
| `stat_off_target_probes.txt` | List of off-target probe IDs |
| `stat_off_target_genes.txt` | List of off-target gene names |
| `stat_missed_probes.txt` | Probes that did not align to their target gene |
| `stat_missed_genes.txt` | Target genes with no aligned probes |
| `track.unmapped.txt` | Probes with no alignments |
| `track.no_hit.txt` | Probes that aligned but passed no acceptance threshold |

---

## Bundled Reference Data

The `data/` directory includes pre-formatted reference files for human:

| Source | Files |
|---|---|
| GENCODE v47 | `data/gencode/gencode.v47.basic.annotation.fmted.fa` / `.gff` |
| RefSeq v110 | `data/refseq/refseq.v110.noAlt.noFix.filtered.fa.gz` / `.gff.gz` |
| CHESS 3.1.3 | `data/chess/chess3.1.3.GRCh38.primary.fa.gz` / `.gff.gz` |

An example gene synonyms file is at `data/gene_synonyms.csv`.

---

## Supported Platforms

- Linux (tested)
- macOS (tested; mummer4 requires Homebrew)

---

## License

See [LICENSE.md](LICENSE.md).
