import argparse
import os
import pyfastx
import pysam
from subprocess import call
import pandas as pd
from datetime import datetime
from enum import Enum
import sys
from Bio.Seq import Seq
import json
from pathlib import Path

# ANSI colors
RED          = '\033[31m'
GREEN        = '\033[32m'
YELLOW       = '\033[33m'
BLUE         = '\033[34m'
CYAN         = '\033[36m'
LIGHT_GREEN  = '\033[92m'
LIGHT_PURPLE = '\033[95m'
RESET        = '\033[0m'

class Mtype(Enum):
    START  = (LIGHT_PURPLE,         "START")
    INFO   = (CYAN,         "PROGRESS")
    RESULT = (GREEN,  "RESULT")
    WARN   = (YELLOW,       "WARNING")
    ERROR  = (RED,          "ERROR")
    DONE   = (LIGHT_PURPLE, "DONE")


# return error message if not in current error list
def message(s, mtype) -> str:
    if mtype not in Mtype:
        raise Exception("Error while printing message")
    return f"{datetime.now()} {mtype.value[0]}{mtype.value[1]}{RESET} {s}"

# check if gff or gtf
def check_annotation_ext(fn) -> str:
    # Strip .gz to get the real extension (e.g. annotation.gff.gz -> .gff)
    name = Path(fn).name.lower()
    if name.endswith('.gz'):
        name = name[:-3]
    file_ext = Path(name).suffix
    att_sep = None
    if file_ext == '.gtf':
        att_sep =  ' '
    elif file_ext in ('.gff', '.gff3'):
        att_sep =  '='
    else:
        print(message(f"unrecognized annotation format '{file_ext}'; must be .gff, .gff3, or .gtf (optionally .gz compressed)", Mtype.ERROR))
        sys.exit(-1)
    return att_sep

# main align function for flipping probes
def align(qfn, tfn, prefix, norc, args) -> str:
    ofn = os.path.join(args.out_dir, f'{prefix}.bam' if args.bam else f'{prefix}.sam')
    if args.binary:
        aligner = args.binary
    else:
        aligner = "bowtie2" if args.bowtie2 else "nucmer"
    print(message(f"running the {aligner} aligner", Mtype.INFO))
    if args.bowtie2: # bt2 flow
        idx_fn = os.path.join(args.out_dir, 'target')
        if not args.skip_index:
            cmd = f'{aligner}-build -q {tfn} {idx_fn} --threads {args.threads}'
            print(cmd)
            call(cmd, shell=True)

        if not os.path.exists(f'{idx_fn}.1.bt2'):
            print(message(f"bt2 index missing; please remove --skip-index flag", Mtype.ERROR))
            sys.exit(-1)

        # add --norc flag if 2nd alignment
        norc_flag = "--norc" if norc else ""

        if args.bam:
            cmd = f'{aligner} -f -a -N 1 --local {norc_flag} -x {idx_fn} ' + \
                f'-U {qfn} --very-sensitive-local --threads {args.threads} ' + \
                f'| samtools view -b -o {ofn} -@ {args.threads}'
        else:
            cmd = f'{aligner} -f -a -N 1 --local {norc_flag} -x {idx_fn} ' + \
                f'-U {qfn} --very-sensitive-local --threads {args.threads} -S {ofn}'
        print(cmd)
        call(cmd, shell=True)
    else: # nucmer flow
        # add -f flag if 2nd alignment
        f_flag = "-f" if norc else ""

        if args.bam:
            temp_sam_fn = os.path.join(args.out_dir, 'temp.sam')
            cmd = f'{aligner} {f_flag} --maxmatch -l {args.min_exact_match} -c 0 -t {args.threads} ' + \
                f'{tfn} {qfn} --sam-long={temp_sam_fn}'
            print(cmd); call(cmd, shell=True)
            cmd = f'samtools view -b -o {ofn} {temp_sam_fn}'
            print(cmd); call(cmd, shell=True)
            cmd = f'rm {temp_sam_fn}'
            print(cmd); call(cmd, shell=True)
        else:
            cmd = f'{aligner} {f_flag} --maxmatch -l {args.min_exact_match} -c 0 -t {args.threads} ' + \
                f'{tfn} {qfn} --sam-long={ofn}'
            print(cmd); call(cmd, shell=True)
    return ofn

def align_nm(qfn, tfn, prefix, args) -> str:
    ofn = os.path.join(args.out_dir, f'{prefix}.mums')
    if args.binary:
        aligner = args.binary
    else:
        aligner = "mummer" # mummer is the only compatible aligner here
    print(message(f"running the {aligner} aligner", Mtype.INFO))
    cmd = f'{aligner} -maxmatch -l {args.min_exact_match} -t {args.threads} ' + \
        f'{tfn} {qfn} > {ofn}'
    print(cmd); call(cmd, shell=True)
    return ofn

def att2dict(s, sep):
    temp = s.split(';')
    d = dict()
    for x in temp:
        kv = x.strip().split(sep)
        if len(kv) != 2: continue
        k = kv[0].strip()
        v = kv[1].strip().replace('"', '')
        d[k] = v
    return d

# tinfo <k,v> = <transcript_id, gene_id>
# build the transcript to gene map
def build_tinfos(fn, att_sep, schema, keep_dot) -> dict:
    df = pd.read_csv(fn, sep='\t', header=None, comment='#')
    df.columns = ['ctg', 'src', 'feat', 'start', 'end', 'score', 'strand', 'frame', 'att']
    # Filter to transcript rows before looping — avoids processing exon/CDS/gene rows
    # (typically 5-10x more rows than transcript rows in a GFF file)
    transcript_df = df[df['feat'] == schema[0]]
    tinfos = dict()
    missing_gene_names = 0
    for row in transcript_df.itertuples(index=False):
        att_d = att2dict(row.att, att_sep)
        if schema[1] not in att_d or schema[2] not in att_d:
            print(message(f"Invalid schema. Expected schema: {schema}. Actual schema: {att_d.keys()}. Change expected schema to correctly locate the necassary information", Mtype.ERROR))
            return None # terminate
        tid = att_d[schema[1]]
        gid = att_d[schema[2]] if keep_dot else att_d[schema[2]].split('.')[0]
        # if you cant get gene name then make None and print out total
        if schema[3] in att_d:
            gname = att_d[schema[3]]
        else:
            gname = None
            missing_gene_names += 1
        if gname:
            temp = gname.split(',')
            if len(temp) > 1:
                temp = [x.strip() for x in temp]
                gname = ';'.join(temp)
        ttype = att_d.get(schema[4], None)
        tinfos[tid] = (gid, gname, ttype)
    print(message(f"missing {missing_gene_names} gene names. If this number is high, the schema may need to be fixed", Mtype.INFO))
    print(message(f"loaded {len(tinfos)} transcripts", Mtype.INFO))

    return tinfos

def write_tinfos(fn, tinfos) -> None:
    with open(fn, 'w') as fh:
        fh.write('transcript_id,gene_id,gene_name,transcript_type\n')
        for x in tinfos:
            y, z, ttype = tinfos[x]
            fh.write(f'{x},{y},{z},{ttype}\n')

def load_tinfos(fn) -> dict:
    df = pd.read_csv(fn)
    # Build dict in one pass using zip — avoids opening the file twice and iterrows overhead
    return dict(zip(df['transcript_id'],
                    zip(df['gene_id'], df['gene_name'], df['transcript_type'])))

def write_lst2file(l, fn) -> None:
    with open(fn, 'w') as fh:
        for x in l:
            fh.write(f'{x}\n')

def read_lst(fn) -> list:
    lst = []
    with open(fn, 'r') as fh:
        for x in fh:
            lst.append(x.strip())
    return lst

def store_params(args, fn):
    with open(fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

def get_unaligned(qfa, ainfos) -> list:
    unaligned = []
    for x in qfa:
        if x.name not in ainfos:
            unaligned.append(x.name)
    return unaligned


# add function to replace gene names with synonyms
def replace_genenames_with_synonyms(query_fasta, gene_syn_csv, out_fasta):
    """
    query_fasta: input FASTA file
    gene_syn_csv: CSV with 2 columns:
                  col1 = original gene name (e.g. WARS),
                  col2 = new gene name (e.g. WARS1)
    out_fasta: output FASTA file with updated headers

    Assumes FASTA headers like:
        >ENSG00000140105|WARS|2532cc5
    and replaces the middle field (WARS) using the synonym map.
    """
    # Read gene synonym mapping
    df = pd.read_csv(gene_syn_csv, header=None)

    # If the CSV has column names, normalize them:
    if df.shape[1] < 2:
        raise ValueError("Synonym file must have at least 2 columns")

    # first 2 columns = original and new
    df = df.iloc[:, :2]

    syn_map = dict(zip(
        df.iloc[:, 0].astype(str).str.strip(),
        df.iloc[:, 1].astype(str).str.strip()
    ))

    # read and update FASTA
    out_lines = []
    replaced_count = 0

    with open(query_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()

                # Split off any whitespace comment after the main ID block
                main_part, *rest_ws = header.split(None, 1)
                rest = rest_ws[0] if rest_ws else ""

                # split using | into 3 parts
                parts = main_part.split('|')

                if len(parts) >= 2:
                    gene_symbol = parts[1]
                    # grab gene name and replace if in syn map
                    if gene_symbol in syn_map:
                        parts[1] = syn_map[gene_symbol]
                        replaced_count += 1

                new_main = "|".join(parts)
                new_header = ">" + new_main
                if rest:
                    new_header += " " + rest
                new_header += "\n"

                out_lines.append(new_header)
            else:
                out_lines.append(line)

    # write new FASTA
    with open(out_fasta, 'w') as f:
        f.writelines(out_lines)
    


