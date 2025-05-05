#!/usr/bin/env python

import pandas as pd
import sys

def parse_brckted_lst(s) -> list:
    trimmed_s = s.strip("[]")
    if not trimmed_s or len(trimmed_s) == 0:
        return None
    brckted_lst = [x.strip() for x in trimmed_s.split(',')]
    return brckted_lst

def load_probes2targets(fn):
    df = pd.read_csv(fn, sep='\t')
    bindings = dict()
    for _, row in df.iterrows():
        pid = row['probe_id']
        n_genes = int(row['n_genes'])
        gene_ids = parse_brckted_lst(row['gene_ids'])
        gene_names = parse_brckted_lst(row['gene_names'])
        cigars = parse_brckted_lst(row['cigars'])
        tx_ids = parse_brckted_lst(row['transcript_ids'])
        tx_types = parse_brckted_lst(row['transcript_types'])
        bindings[pid] = (n_genes, gene_ids, gene_names, cigars, tx_ids, tx_types)
    return bindings

def are_bindings_eq(bndg_1, bndg_2) -> bool:
    if set(bndg_1.keys()) != set(bndg_2.keys()): return False
    for pid in bndg_1:
        b1 = bndg_1[pid]
        b2 = bndg_2[pid]
        if b1[0] != b2[0]: return False
        if set(b1[1]) != set(b2[1]): return False
        if set(b1[2]) != set(b2[2]): return False
        if set(b1[3]) != set(b2[3]): return False
        if set(b1[4]) != set(b2[4]): return False
        if set(b1[5]) != set(b2[5]): return False
    return True

def main(expected, observed, out_fn) -> None:
    b_exp = load_probes2targets(expected)
    b_obs = load_probes2targets(observed)
    is_eq = are_bindings_eq(b_exp, b_obs)
    with open(out_fn, 'w') as fh: fh.write(f'{int(is_eq)}')

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])