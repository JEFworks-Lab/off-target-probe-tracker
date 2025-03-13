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

RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'

class Mtype(Enum):
    PROG = (GREEN, "PROGRESS")
    ERR = (RED, "ERROR")
    WARN = (RED, "WARNING")

def message(s, mtype) -> str:
    if mtype not in Mtype:
        raise Exception("Error while printing message")
    return f"{datetime.now()} {mtype.value[0]}{mtype.value[1]}{RESET} {s}"

def align(qfn, tfn, prefix, norc, args) -> str:
    ofn = os.path.join(args.out_dir, f'{prefix}.bam' if args.bam else f'{prefix}.sam')
    if args.binary:
        aligner = args.binary
    else:
        aligner = "bowtie2" if args.bowtie2 else "nucmer"
    print(message(f"aligner: {aligner}", Mtype.PROG))
    if args.bowtie2: # bt2 flow
        idx_fn = os.path.join(args.out_dir, 'target')
        if not args.skip_index:
            cmd = f'{aligner}-build -q {tfn} {idx_fn} --threads {args.threads}'
            print(cmd)
            call(cmd, shell=True)

        if not os.path.exists(f'{idx_fn}.1.bt2'):
            print(message(f"bt2 index missing; please remove --skip-index flag", Mtype.ERR))
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
    print(message(f"aligner: {aligner}", Mtype.PROG))
    cmd = f'{aligner} -maxmatch -l {args.min_exact_match} -t {args.threads} ' + \
        f'{tfn} {qfn} > {ofn}'
    print(cmd); call(cmd, shell=True)
    return ofn

def att2dict(s, sep):
    temp = s.split(';')
    d = dict()
    for x in temp:
        kv = x.split(sep)
        if len(kv) != 2: continue
        k = kv[0].strip()
        v = kv[1].strip()
        d[k] = v
    return d

# tinfo <k,v> = <transcript_id, gene_id>
def build_tinfos(fn, att_sep, schema, keep_dot) -> dict:
    df = pd.read_csv(fn, sep='\t', header=None)
    df.columns = ['ctg', 'src', 'feat', 'start', 'end', 'score', 'strand', 'frame', 'att']
    tinfos = dict()
    ctr = 0
    for _, row in df.iterrows():
        if row['feat'] == schema[0]:
            ctr += 1
            att_d = att2dict(row['att'], att_sep)
            if schema[1] not in att_d or schema[2] not in att_d: 
                print(message(f"Invalid schema", Mtype.ERR))
                return None # terminate
            tid = att_d[schema[1]]
            gid = att_d[schema[2]] if keep_dot else att_d[schema[2]].split('.')[0]
            gname = att_d[schema[3]] if schema[3] in att_d else None
            if gname:
                temp = gname.split(',')
                if len(temp) > 1:
                    temp = [x.strip() for x in temp]
                    gname = ';'.join(temp)
            ttype = att_d.get(schema[4], None) # Caleb: add transcript type
            tinfos[tid] = (gid, gname, ttype) # Caleb: add transcript type
    print(message(f"loaded {ctr} transcripts", Mtype.PROG))
    return tinfos

def write_tinfos(fn, tinfos) -> None:
    with open(fn, 'w') as fh:
        fh.write('transcript_id,gene_id,gene_name,transcript_type\n') # Caleb: add transcript type
        for x in tinfos:
            y, z, ttype = tinfos[x] # Caleb: add transcript type
            fh.write(f'{x},{y},{z},{ttype}\n') # Caleb: add transcript type

def load_tinfos(fn) -> dict:
    tinfos = dict()
    df = pd.read_csv(fn)
    with open(fn, 'r') as fh:
        for _, row in df.iterrows():
            tinfos[row['transcript_id']] = (row['gene_id'], row['gene_name'], row['transcript_type']) # Caleb: add transcript type
    return tinfos

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