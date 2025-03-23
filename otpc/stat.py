from otpc.commons import *
import pandas as pd
from collections import namedtuple

class Hit:
    def __init__(self, gid, gname, tid, cigar, ttype):
        self.gid = gid
        self.gname = gname
        self.tid = tid
        self.cigar = cigar
        self.ttype = ttype

class Probe:
    def __init__(self, id: str, gids: list, gnames: list, \
                n_hits: int, tids: list, cigars: list, ttypes: list):
        self.id = id
        self.n_hits = n_hits
        self.hits = [Hit(gid, gname, tid, cigar, ttype) for \
                    gid, gname, tid, cigar, ttype in \
                    zip(gids, gnames, tids, cigars, ttypes)]
        assert self.n_hits == len(self.hits) # sanity check

class ProbeGene:
    def __init__(self, id: str, name: str):
        self.id = id # target gene ID
        self.name = name # target gene name
        self.n = 0 # number of probes
        self.probe_hits = dict()
    
    def add_prb(self, p):
        self.probe_hits[p.id] = p
        self.n += 1
    
def parse_brckted_lst(s) -> list:
    trimmed_s = s.strip("[]")
    if not trimmed_s or len(trimmed_s) == 0:
        return None
    brckted_lst = [x.strip() for x in trimmed_s.split(',')]
    return brckted_lst

# NOTE: custom to 10x probe IDs
def get_probe_info(s) -> tuple:
    temp = s.split("|")
    if len(temp) != 3: 
        return None
    p_gid, p_gname, pid = temp
    return (p_gid, p_gname, pid)

def is_all_ps(l) -> bool:
    return all('pseudogene' in x for x in l)

def is_all_non_pc(l) -> bool:
    return all(x != 'protein_coding' and x != 'mRNA' for x in l)

def load_track_results(df, d, gene_syns, exclude_ps, pc_only) -> dict:
    prb_gene_tbl = dict()
    missed_target_prbs = []
    multi_target_prbs = []
    for i, row in df.iterrows():
        pinfo = get_probe_info(row['probe_id'])
        if not pinfo:
            print(message(f"failed to load probe info from line #{i}", Mtype.WARN))
            continue
        p_gid, p_gname, pid = pinfo
        cigars = parse_brckted_lst(row['cigars'])
        ttypes = parse_brckted_lst(row['transcript_types'])
        tids = parse_brckted_lst(row['transcript_ids'])
        gnames = parse_brckted_lst(row['gene_names'])
        gids = parse_brckted_lst(row['gene_ids'])
        assert len(cigars) == len(ttypes) == len(tids) == len(gnames) == len(gids) # sanity check
        prb = Probe(pid, gids, gnames, len(cigars), tids, cigars, ttypes)

        # count missed_target and off_target probes
        off = False
        missed = True
        p_genes = [p_gname]
        if p_gname in gene_syns:
            p_genes.append(gene_syns[p_gname])
        for i, x in enumerate(gnames):
            if exclude_ps and 'pseudogene' in ttypes[i]: continue
            if pc_only and (ttypes[i] != 'protein_coding' and ttypes[i] != 'mRNA'): continue

            if x in p_genes:
                missed = False
            else:
                off = True
        if off:
            multi_target_prbs.append(row['probe_id'])
        if missed:
            missed_target_prbs.append(row['probe_id'])

        # CAUTION: there might be >1 distinct gene_ids for the same gene_name
        if p_gname not in prb_gene_tbl:
            prb_gene_tbl[p_gname] = ProbeGene(p_gid, p_gname)
        prb_gene_tbl[p_gname].add_prb(prb)
        assert pid in prb_gene_tbl[p_gname].probe_hits # sanity check
    print(message(f"number of probes missing targets: {len(missed_target_prbs)}", Mtype.PROG))
    write_lst2file(missed_target_prbs, os.path.join(d, 'stat_missed_probes.txt'))
    print(message(f"number of probes with off-target binding: {len(multi_target_prbs)}", Mtype.PROG))
    write_lst2file(multi_target_prbs, os.path.join(d, 'stat_off_target_probes.txt'))
    return prb_gene_tbl

def summarize(prb_gene_tbl, exclude_ps, pc_only) -> dict:
    agg = dict()
    for p_gname in prb_gene_tbl:
        temp = dict()
        for pid in prb_gene_tbl[p_gname].probe_hits:
            prb = prb_gene_tbl[p_gname].probe_hits[pid]
            temp2 = set()
            for ht in prb.hits:
                if 'pseudogene' in ht.ttype and exclude_ps:
                    continue
                # TODO: allow user to specify pattern for pc txes
                if pc_only and (ht.ttype != 'protein_coding' and ht.ttype != 'mRNA'):
                    continue
                if ht.gname not in temp:
                    temp[ht.gname] = [1, 1, {ht.gid}]
                    temp2.add(ht.gname)
                else:
                    if ht.gname not in temp2:
                        temp[ht.gname][1] += 1
                        temp2.add(ht.gname)
                    temp[ht.gname][0] += 1
                    temp[ht.gname][2].add(ht.gid)
        agg[p_gname] = temp
    return agg

def write_summary(d, agg, pgene_info, gene_syns):
    fn = os.path.join(d, 'stat.summary.tsv')
    clpsed = dict()
    with open(fn, 'w') as fh:
        fh.write('target_gene\tn\tgene_name\tn_hits\tn_probes\tgene_ids\n')
        for p_gname in agg:
            clpsed[p_gname] = [[], [], []]
            for gname in agg[p_gname]:
                temp = agg[p_gname][gname]
                clpsed[p_gname][0].append(gname)
                clpsed[p_gname][1].append(str(temp[0]))
                clpsed[p_gname][2].append(str(temp[1]))
                fh.write(f'{p_gname}\t{pgene_info[p_gname]}\t{gname}')
                fh.write(f'\t{temp[0]}\t{temp[1]}\t[{','.join(temp[2])}]\n')

    missed_target = []
    off_target = []
    fn = os.path.join(d, 'collapsed_summary.tsv')
    print(message(f"{len(clpsed)} / {len(pgene_info)} probe genes with at least 1 probe binding", Mtype.PROG))
    with open(fn, 'w') as fh:
        fh.write('target_gene\tn\taligned_to\tn_hits\tn_probes\n')
        for p_gname in clpsed:
            temp = clpsed[p_gname]
            fh.write(f'{p_gname}\t{pgene_info[p_gname]}\t[{','.join(temp[0])}]')
            fh.write(f'\t[{','.join(temp[1])}]\t[{','.join(temp[2])}]\n')
            gnames = temp[0]

            off = False
            missed = True
            p_genes = [p_gname]
    
            if p_gname in gene_syns:
                p_genes.append(gene_syns[p_gname])
            
            for x in gnames:
                if x in p_genes:
                    missed = False
                else:
                    off = True

            if off:
                off_target.append(p_gname)
            if missed:
                missed_target.append(p_gname)
    print(message(f"number of missed probe genes: {len(missed_target)}", Mtype.PROG))
    print(message(f"number of off-target probe genes: {len(off_target)}", Mtype.PROG))
    write_lst2file(missed_target, os.path.join(d, 'stat_missed_genes.txt'))
    write_lst2file(off_target, os.path.join(d, 'stat_off_target_genes.txt'))

def load_pgene_info(fn) -> dict:
    pgene_info = dict()
    fa = pyfastx.Fasta(fn)
    for x in fa:
        pinfo = get_probe_info(x.name)
        if not pinfo:
            print(message(f"failed to load probe info", Mtype.WARN))
            continue
        _, p_gname, _ = pinfo
        if p_gname not in pgene_info:
            pgene_info[p_gname] = 1
        else:
            pgene_info[p_gname] += 1
    return pgene_info

# TODO: is this optimal?
def load_gene_syns(fn) -> dict:
    gene_syns = dict()
    with open(fn, 'r') as fh:
        for ln in fh:
            pair = ln.strip().split(',')
            gene_syns[pair[0]] = pair[1]
    return gene_syns

def main(args) -> None:
    gene_syns = load_gene_syns(args.syn_file) if args.syn_file else []
    pgene_info = load_pgene_info(args.query)
    track_df = pd.read_csv(args.in_file, sep='\t')
    prb_gene_tbl = load_track_results(track_df, args.out_dir, gene_syns, args.exclude_pseudo, args.pc_only)
    print(message(f"number of probe genes: {len(prb_gene_tbl)}", Mtype.PROG))
    agg = summarize(prb_gene_tbl, args.exclude_pseudo, args.pc_only)
    write_summary(args.out_dir, agg, pgene_info, gene_syns)