from opt.commons import *

# check for missing gene IDs and names
def check_sinfo_completeness(qinfos, tinfos):
    tgt_gids = set([x[0] for x in tinfos.values()])
    tgt_gnames = set([x[1] for x in tinfos.values()])
    missing_gids = set()
    missing_gnames = set()
    for qname in qinfos:
        gid, gname = qinfos[qname]
        if gid not in tgt_gids:
            missing_gids.add(gid)
        if gname not in tgt_gnames:
            missing_gnames.add(gname)
    print(message(f"{len(missing_gids)} target gene IDs missing from the source. If reference genome is not GENCODE, ignore this. See missing_in_src_ids.txt", Mtype.RESULT))
    print(message(f"{len(missing_gnames)} target gene names missing from the source. See missing_in_src_names.txt ", Mtype.RESULT))
    return missing_gids, missing_gnames

# read in probe file
def load_pinfos(fn):
    qfa = pyfastx.Fasta(fn)
    qinfos = dict()
    for q in qfa:
        temp = q.name.split("|")
        gid = temp[0]
        gname = temp[1]
        qinfos[q.name] = (gid, gname)
    return qinfos, qfa

def write_pinfos(out_dir, qinfos) -> None:
    fn = os.path.join(out_dir, 'probe_infos.csv')
    with open(fn, 'w') as fh:
        fh.write('probe_id,gene_id,gene_name\n')
        for x in qinfos:
            fh.write(f'{x},{qinfos[x][0]},{qinfos[x][1]}\n')

# read in sam file
def load_bam(fn, qinfos, tinfos, is_bam) -> dict:
    ainfos = dict()
    fh = pysam.AlignmentFile(fn, 'rb') if is_bam else pysam.AlignmentFile(fn, 'r')
    for brec in fh:
        qname = brec.query_name # qname example = 'ENSG00000170458|CD14|22f9405'
        if qname not in ainfos: # ainfos example = 'ENSG00000170458|CD14|22f9405': [False, False, False, False]}
            ainfos[qname] = []
        if brec.is_unmapped: continue
        tname = brec.reference_name # get transcript name
        assert qname in qinfos # sanity check
        target_gid, target_gname = qinfos[qname]
        # if target_gid == tinfos[tname][0]: # using gene name instead of gene id since refseq/chess dont have genecode gene id
        if target_gname == tinfos[tname][1]:
            ainfos[qname].append(brec.is_forward)
    fh.close()
    return ainfos # dict of query name and list of T/F if gene is forward

def flip(ainfos, pfa, out_dir):
    missing_origin = []
    # Use a set for to_rc so membership checks during FASTA writing are O(1)
    # instead of O(n) — critical when there are thousands of probes
    to_rc = set()
    for qname in ainfos:
        if len(ainfos[qname]) == 0:
            missing_origin.append(qname)
        elif all(ainfos[qname]):
            continue
        else:
            assert all(not x for x in ainfos[qname])
            to_rc.add(qname)
    print(message(f"{len(missing_origin)} / {len(pfa)} probes not mapped to their origin. See flip.missing_ori.txt", Mtype.RESULT))
    print(message(f"{len(to_rc)} / {len(pfa)} probes to flip (i.e., reverse complement). See rev_cmped_probes.txt", Mtype.RESULT))
    fn = os.path.join(out_dir, 'fwd_oriented.fa')
    with open(fn, 'w') as fh:
        for q in pfa:
            if q.name in to_rc:
                seq = Seq(q.seq)
                out_s = seq.reverse_complement()
            else:
                out_s = q.seq
            fh.write(f'>{q.name}\n{out_s}\n')
    return missing_origin, list(to_rc)

def main(args) -> None:
    print(message(f"FLIP module is aligning input probes to source transcripts", Mtype.START))
    bfn = align(args.query, args.target, 'flip', False, args)
    # att_sep = ' ' if args.gtf else '='
    # check if gtf or gff
    att_sep = check_annotation_ext(args.annotation)
    # save patch for flip genes
    fn = os.path.join(args.out_dir, 'flip_t2g.csv')
    print(message(f"loading source transcriptome infos", Mtype.INFO))
    # get transcript to gene mapping
    if not os.path.exists(fn) or args.force:
        print(message(f"building transcript to gene mappings", Mtype.INFO))
        sinfos = build_tinfos(args.annotation, att_sep, args.schema, args.keep_dot) # sinfos example = dict of 'ENST00000382296.4': ('ENSG00000205916', 'DAZ4', 'protein_coding')
        write_tinfos(fn, sinfos)
    else:
        sinfos = load_tinfos(fn)
    print(message(f"loading input probes infos", Mtype.INFO))
    pinfos, pfa = load_pinfos(args.query) # pinfos example = dict of 'ENSG00000285437|POLR2J3|b33269d': ('ENSG00000285437', 'POLR2J3')
    
    # save missing ids and genes
    missing_gids, missing_gnames = check_sinfo_completeness(pinfos, sinfos)
    if len(missing_gids) > 0:
        write_lst2file(missing_gids, os.path.join(args.out_dir, 'missing_in_src_ids.txt'))
    if len(missing_gnames) > 0:
        write_lst2file(missing_gnames, os.path.join(args.out_dir, 'missing_in_src_names.txt'))
    
    # write flip files
    print(message(f"parsing alignment results", Mtype.INFO))
    ainfos = load_bam(bfn, pinfos, sinfos, args.bam)
    unaligned = get_unaligned(pfa, ainfos)
    print(message(f"{len(unaligned)} / {len(pfa)} probes unmapped. See flip.unmapped.txt", Mtype.RESULT))
    missing_origin, flipped = flip(ainfos, pfa, args.out_dir)
    if len(unaligned) > 0:
        write_lst2file(unaligned, os.path.join(args.out_dir, 'flip.unmapped.txt'))
    if len(missing_origin) > 0:
        write_lst2file(missing_origin, os.path.join(args.out_dir, 'flip.missing_ori.txt'))    
    if len(flipped) > 0:
        write_lst2file(flipped, os.path.join(args.out_dir, 'rev_cmped_probes.txt'))
    print(message(f"wrote forward oriented probes to a file", Mtype.INFO))
    print(message(f"finished flipping probes", Mtype.DONE))