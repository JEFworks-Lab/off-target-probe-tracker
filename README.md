# OPT — Off-target Probe Tracker

OPT identifies potential off-target binding of probe sequences against a reference transcriptome using nucleotide alignment (nucmer). The goal of OPT is to help evaluate probe specificity before experiments by detecting probes that may hybridize to unintended transcripts.


Hallinan et al., *eLife* 2025. https://elifesciences.org/reviewed-preprints/107070


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
4. Decompress examples probes and all reference annotations

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
2. **Input Files** — provide paths to your probe FASTA and reference files. Use the Browse buttons or type paths directly.
   - Select an **annotation format** preset: **GENCODE**, **RefSeq**, **CHESS**, or **Other** (custom schema). The correct GFF/GTF schema is applied automatically.
   - Select **All (GENCODE + CHESS + RefSeq)** to run against all three reference annotations sequentially and merge the results into a unified off-target table.
   - An optional **gene synonyms CSV** can be provided to remap gene names that differ between your probe FASTA and the reference (e.g. `WARS` → `WARS1`).
3. **Analysis Options:**
   - **Pad length** — number of bases at each probe end where mismatches are tolerated (default: off). Used in the original paper for Xenium probes, which are circular and can tolerate terminal mismatches.
   - **Max mismatches anywhere** — allow up to N mismatches anywhere in the full probe sequence (default: off). Can be combined with pad length: when both are set, both conditions must be satisfied.
4. Click **Run OPT** to run all three modules (flip → track → stat) and view results in the dashboard below.

The results dashboard shows:
- **Brief Summary** - total genes with off-target binding, total genes with protein-coding off-targets, total probes with off-target binding
- **Gene-level off-target table** — one row per target gene → off-target gene pair, with biotype badges, CIGAR strings, and source annotation. Filterable by biotype and sortable by any column.
- **Probe-level detail table** (expandable) — one row per probe, showing off-target genes, biotypes, and CIGAR strings (consistent counts, `|`-delimited).
- **Download buttons** for all key output files.

To load results from a previous run without re-running OPT, set the Output Directory to your previous run folder and click **Load previous results**.

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
>gene_id|gene_name|unique_id
```

Example:
```
>ENSG00000170458|CD14|22f9405
ATCGATCGATCGATCGATCG...
```

### Target transcript FASTA

Standard nucleotide FASTA of transcript sequences (.fa or .fasta). We recommend extracting these with [gffread](https://github.com/gpertea/gffread):

```bash
gffread -w transcripts.fa -g genome.fa annotation.gff
```

> **Note:** The web app requires uncompressed `.fa` or `.fasta` files. The CLI accepts any format that nucmer/Bowtie2 supports.

### Annotation GFF/GTF

Standard GFF3 or GTF format (.gff, .gff3, or .gtf). GENCODE, RefSeq, and CHESS formats are all supported. Select the matching preset in the web app, or use `--schema` on the command line for non-standard formats.

> **Note:** The web app requires uncompressed annotation files. Gzip-compressed files (`.gz`) are supported via the CLI only.

### Gene Synonyms CSV (optional)

Two-column CSV mapping probe gene names to annotation gene names. No header required:

```
WARS,WARS1
CARS,CARS1
```

Use this when gene names in your probe FASTA differ from those in the reference annotation.

---

## GFF/GTF Schema

The `--schema` argument specifies five comma-separated field names used to parse the annotation. Built-in presets for common formats:

| Format | Schema string |
|---|---|
| GENCODE GFF | `transcript,ID,Parent,gene_name,transcript_type` |
| RefSeq GFF | `transcript,ID,Parent,gene,gbkey` |
| CHESS GFF | `transcript,ID,Parent,gene_name,gene_type` |
| GTF (general) | `transcript,transcript_id,gene_id,gene_name,transcript_type` |

| Position | Description |
|---|---|
| 1 | Feature type (column 3 of the GFF/GTF) |
| 2 | Transcript ID attribute |
| 3 | Parent gene attribute |
| 4 | Gene name attribute |
| 5 | Transcript type / biotype attribute |

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

When running in **All (GENCODE + CHESS + RefSeq)** mode, each annotation runs in its own subdirectory (`gencode/`, `chess/`, `refseq/`) and results are merged into the base output directory with an added `reference_annotation` column.

---

## Bundled Reference Data

The `data/` directory includes pre-formatted reference files for human (GRCh38):

| Source | Files |
|---|---|
| GENCODE v47 | `data/gencode/gencode.v47.basic.annotation.fmted.fa` / `.gff` |
| RefSeq v110 | `data/refseq/refseq.v110.noAlt.noFix.filtered.fa.gz` / `.gff.gz` |
| CHESS 3.1.3 | `data/chess/chess3.1.3.GRCh38.primary.fa.gz` / `.gff.gz` |

An example gene synonyms file is at `data/gene_synonyms.csv`.

The `install.sh` script automatically decompresses all `.gz` files in the `data/` directory. To decompress manually:

```bash
find data/ -name "*.gz" -exec gunzip -k {} \;
```

---

## Supported Platforms

- Linux (tested)
- macOS (tested; mummer4 requires Homebrew)

---

## License

See [LICENSE.md](LICENSE.md).
