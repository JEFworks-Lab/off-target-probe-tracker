from opt.commons import *

def convert_md2bit(s):
    # Use list + join instead of repeated string concatenation (O(n) vs O(n²))
    parts = []
    running = ""
    for c in s:
        if c.isdigit():
            running += c
        else:
            if running:
                parts.append('1' * int(running))
            parts.append('0')
            running = ""
    if running:
        parts.append('1' * int(running))
    return ''.join(parts)

def convert_md2bit_nucmer(s, tstart):
    parts = []
    running = ""
    for c in s:
        if c.isdigit():
            running += c
        else:
            if running:
                parts.append('1' * int(running))
            parts.append('0')
            running = ""
    if running:
        parts.append('1' * (int(running) - tstart))
    return ''.join(parts)

def convert_md2bit_del(s):
    parts = []
    pos = 0  # tracks current output length to record mismatch positions
    running = ""
    ignore = False
    mismatch_info = []
    for c in s:
        if c.isdigit():
            ignore = False
            running += c
        else:
            if ignore: continue
            if c == '^':
                if running:
                    n = int(running)
                    parts.append('1' * n)
                    pos += n
                running = ""
                ignore = True
                continue
            else:
                if running:
                    n = int(running)
                    parts.append('1' * n)
                    pos += n
                mismatch_info.append(pos)
                parts.append('0')
                pos += 1
                running = ""
    if running:
        parts.append('1' * int(running))
    return ''.join(parts), mismatch_info

def convert_md2bit_nucmer_del(s, tstart):
    parts = []
    pos = 0
    running = ""
    ignore = False
    mismatch_info = []
    for c in s:
        if c.isdigit():
            ignore = False
            running += c
        else:
            if ignore: continue
            if c == '^':
                if running:
                    n = int(running)
                    parts.append('1' * n)
                    pos += n
                running = ""
                ignore = True
                continue
            else:
                if running:
                    n = int(running)
                    parts.append('1' * n)
                    pos += n
                mismatch_info.append(pos)
                parts.append('0')
                pos += 1
                running = ""
    if running:
        parts.append('1' * (int(running) - tstart))
    return ''.join(parts), mismatch_info

def convert_cigar2bit(tup):
    bit_s = ""
    left_clip = 0
    right_clip = 0
    ins_info = []
    for x in tup:
        op, l = x
        if op == 0: # match
            bit_s += '1' * l
        elif op == 1 or op == 4: # soft clip or ins
            if op == 4:
                if len(bit_s) == 0:
                    left_clip = l
                else:
                    right_clip = l
            if op == 1:
                ins_info.append((bit_s.count('1'), l))
            bit_s += '0' * l
    return bit_s, (left_clip, right_clip), ins_info

def convert_cigar2bit_del(tup, n, mismatch_info):
    bit_s_lst = []
    bit_s = ""
    running = 0
    mut_ctr = 0
    for x in tup:
        op, l = x
        if op == 0: # match
            temp = '1' * l
            # TODO: test if this behaves as expected
            # switch '1' to '0' if a mismatch is reported in this stretch of matches
            if len(mismatch_info) > 0 and mut_ctr != len(mismatch_info):
                temp_lst = list(temp)
                for m in mismatch_info:
                    if m >= running and m < running + l:
                        temp_lst[m - running] = '1'
                        mut_ctr += 1
                temp = ''.join(temp_lst)
            running += l
            bit_s += temp
        elif op == 1 or op == 4: # soft clip or ins
            bit_s += '0' * l
        elif op == 2:
            temp = '0' * len(bit_s)
            bit_s += '0' * (n - len(bit_s))
            bit_s_lst.append(bit_s)
            bit_s = temp
    bit_s_lst.append(bit_s)
    return bit_s_lst

def bitwise_and(s1, s2):
    assert len(s1) == len(s2)
    return ''.join('1' if a == '1' and b == '1' else '0' for a, b in zip(s1, s2))

def char2sym(char):
    if char == '0':
        return 'X'
    return '='

def compress_bvec(bvec):
    out = []
    curr_char = bvec[0]
    ctr = 1
    for char in bvec[1:]:
        if char == curr_char:
            ctr += 1
        else:
            out.append(f"{char2sym(curr_char)}{ctr}")
            curr_char = char
            ctr = 1
    out.append(f"{char2sym(curr_char)}{ctr}")
    return ''.join(out)

def track_target_pad(fn, qfa, pad, tinfos, is_nucmer, max_mismatches=-1) -> dict:
    ainfos = dict()
    # Cache per-probe length and critical bit vector — each probe can appear in thousands
    # of alignment records, so recomputing these every time is wasteful
    qlen_cache = {}
    crit_cache = {}
    with pysam.AlignmentFile(fn, 'rb') as fh:
        for brec in fh:
            qname = brec.query_name
            # NOTE: nucmer doesn't output unmapped reads
            if brec.is_unmapped:
                continue
            elif brec.is_supplementary:
                continue
            else:
                tname = brec.reference_name
                if qname not in qlen_cache:
                    qlen_cache[qname] = len(qfa[qname].seq)
                    crit_bvec = "0" * pad + "1" * (qlen_cache[qname] - 2 * pad) + "0" * pad
                    assert len(crit_bvec) == qlen_cache[qname] # sanity check
                    crit_cache[qname] = int(crit_bvec, 2)
                qlen = qlen_cache[qname]
                crit_dvec = crit_cache[qname]
                if qname not in ainfos:
                    ainfos[qname] = set() # empty if no brec is passing
                    # NOTE: this used to be a set; explains the discrepancy in the output
                cigar = brec.cigarstring
                cigar_tups = brec.cigartuples
                num_mismatch = int(brec.get_tag('NM'))
                md_tag = brec.get_tag('MD')
                tstart = brec.reference_start
                if cigar == f'{qlen}M':
                    if num_mismatch == 0:
                        ainfos[qname].add((tname, (tinfos[tname][0], tinfos[tname][1]), \
                                        tinfos[tname][2], f'={qlen}'))
                    else:
                        if is_nucmer:
                            md_bvec = convert_md2bit_nucmer(md_tag, tstart)
                        else:
                            md_bvec = convert_md2bit(md_tag)
                        crit_pass = crit_dvec & int(md_bvec, 2) == crit_dvec
                        if max_mismatches >= 0:
                            accept = num_mismatch <= max_mismatches and (pad == 0 or crit_pass)
                        else:
                            accept = crit_pass
                        if accept:
                            ainfos[qname].add((tname, (tinfos[tname][0], tinfos[tname][1]), \
                                        tinfos[tname][2], compress_bvec(md_bvec)))
                else:
                    if 'D' in cigar: # handle deletions separately
                        # NOTE: no need to check num_mismatch == 0 as dels count as mismatches (i.e., nm > 0 guaranteed)
                        if is_nucmer:
                            md_bvec, mismatch_info = convert_md2bit_nucmer_del(md_tag, tstart)
                        else:
                            md_bvec, mismatch_info = convert_md2bit_del(md_tag)
                        cigar_bvecs = convert_cigar2bit_del(cigar_tups, qlen, mismatch_info)
                        hit = False
                        final_bvec = None
                        for bvec in cigar_bvecs:
                            if crit_dvec & int(bvec, 2) == crit_dvec:
                                final_bvec = bvec
                                hit = True
                                break
                        if max_mismatches >= 0:
                            accept = num_mismatch <= max_mismatches and (pad == 0 or hit)
                        else:
                            accept = hit
                        if accept:
                            if final_bvec is None:
                                final_bvec = cigar_bvecs[0] if cigar_bvecs else '1' * qlen
                            ainfos[qname].add((tname, (tinfos[tname][0], tinfos[tname][1]), \
                                        tinfos[tname][2], compress_bvec(final_bvec)))
                    else:
                        cigar_bvec, clip_info, ins_info = convert_cigar2bit(cigar_tups)
                        if num_mismatch == 0: # accounts for cases with just soft clips
                            final_bvec = cigar_bvec # NOTE: ins counts as a mismatch
                        else:
                            if is_nucmer:
                                md_bvec = convert_md2bit_nucmer(md_tag, tstart)
                            else:
                                md_bvec = convert_md2bit(md_tag)
                            if len(ins_info) > 0:
                                # convert md_bvec to account for soft-clipped and inserted bases
                                temp = ""
                                prev_p = None
                                for i in range(len(ins_info)):
                                    p, l = ins_info[i]
                                    if i == 0:
                                        temp += md_bvec[:p] + '0' * l
                                    else:
                                        temp += md_bvec[prev_p:p] + '0' * l
                                    prev_p = p
                                temp += md_bvec[ins_info[-1][0]:]
                                temp = '0' * clip_info[0] + temp + '0' * clip_info[1]
                                assert len(temp) == len(cigar_bvec) # sanity check
                                final_bvec = bitwise_and(cigar_bvec, temp)
                            else:
                                temp = '0' * clip_info[0] + md_bvec + '0' * clip_info[1]
                                assert len(temp) == len(cigar_bvec) # sanity check
                                final_bvec = bitwise_and(cigar_bvec, temp)
                        crit_pass = crit_dvec & int(final_bvec, 2) == crit_dvec
                        if max_mismatches >= 0:
                            # Use bit-vector 0-count: NM tag excludes soft-clips/insertions
                            effective_nm = final_bvec.count('0')
                            accept = effective_nm <= max_mismatches and (pad == 0 or crit_pass)
                        else:
                            accept = crit_pass
                        if accept:
                            ainfos[qname].add((tname, (tinfos[tname][0], tinfos[tname][1]), \
                                            tinfos[tname][2], compress_bvec(final_bvec)))
    return ainfos

def load_mums(fn) -> dict:
    mums = dict()
    with open(fn, 'r') as fh:
        for ln in fh:
            clean_ln = ln.strip()
            if clean_ln[0] == '>':
                qname = clean_ln.split()[1].replace("> ", "")
            else:
                temp = clean_ln.split()
                if qname not in mums:
                    mums[qname] = [(temp[0], int(temp[1]), int(temp[2]), int(temp[3]))]
                else:
                    mums[qname].append((temp[0], int(temp[1]), int(temp[2]), int(temp[3])))
    return mums

def check_lft_and_rgt(mrec, qseq, qlen, tgt_fa, max_nm):
    # Accepts pre-fetched qseq and qlen — caller caches these per probe to avoid
    # repeated pyfastx index lookups when the same probe has multiple alignment hits
    tname, tst, qst, mlen = mrec
    if mlen == 40:
        return (True, '1' * mlen, 0)
    qst -= 1
    tst -= 1
    qen = qst + mlen
    ten = tst + mlen
    lft_qos = qst
    rgt_qos = qlen - qen
    tseq = tgt_fa[tname].seq
    lft_tseq = tseq[tst - lft_qos:tst]
    rgt_tseq = tseq[ten:ten + rgt_qos]
    lft_qseq = qseq[qst - lft_qos:qst]
    rgt_qseq = qseq[qen:qen + rgt_qos]
    if len(lft_tseq) != lft_qos: # tseq runs out at 5'
        return (False, None, -1)
    if len(rgt_tseq) != rgt_qos: # tseq runs out at 3'
        return (False, None, -1)
    # Build match vectors with list comprehension + join instead of per-char concatenation
    lft_bits = ['1' if lft_tseq[i] == lft_qseq[i] else '0' for i in range(lft_qos)]
    lft_nm = lft_bits.count('0')
    rgt_bits = ['1' if rgt_tseq[i] == rgt_qseq[i] else '0' for i in range(rgt_qos)]
    rgt_nm = rgt_bits.count('0')
    mvec = ''.join(lft_bits) + ('1' * mlen) + ''.join(rgt_bits)
    assert len(mvec) == qlen # sanity check
    nm = lft_nm + rgt_nm
    return (nm <= max_nm, mvec, nm)

def track_target_one_mismatch(fn, qfa, tfa, tinfos) -> dict:
    mums = load_mums(fn)
    ainfos = dict()
    for qname, mrecs in mums.items():
        ainfos[qname] = set()
        # Pre-fetch probe sequence once per probe — each probe can appear in many
        # alignment records (mrecs), so caching avoids repeated pyfastx index hits
        q = qfa[qname]
        qseq = q.seq
        qlen = len(qseq)
        for mrec in mrecs:
            tname = mrec[0]
            is_pass, mvec, _ = check_lft_and_rgt(mrec, qseq, qlen, tfa, 1)
            if is_pass:
                ainfos[qname].add((tname, (tinfos[tname][0], tinfos[tname][1]), \
                                            tinfos[tname][2], compress_bvec(mvec)))
    return ainfos
    
def _summarize_types(types):
    """Return the most informative type label for a set of transcript types for one gene."""
    if 'protein_coding' in types or 'mRNA' in types:
        return 'protein_coding'
    return ';'.join(sorted(types))


def write_results(ainfos, d) -> list:
    # Single pass over ainfos — build both output files simultaneously to avoid
    # iterating the dict twice and extracting fields twice per probe
    header    = 'probe_id\tn_genes\tgene_ids\tgene_names\tcigars\ttranscript_ids\ttranscript_types\n'
    header_ot = ('probe_id\tn_genes\tgene_ids\tgene_names\tcigars\ttranscript_ids\ttranscript_types'
                 '\tprobe_gene\tofftarget_gene_names\tofftarget_gene_types\tconcern_level\n')
    fn    = os.path.join(d, 'probe2targets.tsv')
    fn_ot = os.path.join(d, 'probe2targets_offtargets.tsv')
    no_hit = []
    with open(fn, 'w') as fh, open(fn_ot, 'w') as fh_ot:
        fh.write(header)
        fh_ot.write(header_ot)
        for qname, hits in ainfos.items():
            if len(hits) == 0:
                no_hit.append(qname)
                continue
            tnames = [x[0] for x in hits]
            genes  = [x[1] for x in hits]
            gids   = [x[0] for x in genes]
            gnames = [x[1] for x in genes]
            ttypes = [x[2] for x in hits]
            cigars = [x[3] for x in hits]
            try:
                assert len(gids) == len(gnames) # sanity check
            except:
                print(message(f">1 reference gene IDs might share the same gene name", Mtype.WARN))
                print(gids)
                print(gnames)

            gids_s   = ','.join(str(x) for x in gids)
            gnames_s = ','.join('None' if x is None else str(x) for x in gnames)
            cigar_s  = ','.join(str(x) for x in cigars)
            ttypes_s = ','.join(str(x) for x in ttypes)
            tnames_s = ','.join(str(x) for x in tnames)
            n_genes  = len(set(gnames))
            row = f'{qname}\t{n_genes}\t[{gids_s}]\t[{gnames_s}]\t[{cigar_s}]\t[{tnames_s}]\t[{ttypes_s}]\n'
            fh.write(row)
            if n_genes > 1:
                # Derive probe's intended gene from its ID (format: gid|gname|accession)
                parts = qname.split('|')
                probe_gene = parts[1] if len(parts) >= 2 else qname

                # Build gene_name → set of transcript types
                gene_to_types = {}
                for gname_i, ttype_i in zip(gnames, ttypes):
                    gene_to_types.setdefault(gname_i, set()).add(str(ttype_i))

                # Unique off-target gene names (not the probe's target), order preserved
                seen_ot = set()
                ot_genes = []
                for gname_i in gnames:
                    if gname_i != probe_gene and gname_i not in seen_ot:
                        seen_ot.add(gname_i)
                        ot_genes.append(gname_i)

                ot_gene_names_s = ','.join(str(g) for g in ot_genes)
                ot_gene_types_s = ','.join(_summarize_types(gene_to_types[g]) for g in ot_genes)

                # Concern level based on off-target gene types
                ot_type_sets = [gene_to_types[g] for g in ot_genes]
                if any('protein_coding' in ts or 'mRNA' in ts for ts in ot_type_sets):
                    concern = 'high'
                elif any(not any('pseudogene' in t for t in ts) for ts in ot_type_sets):
                    concern = 'medium'
                else:
                    concern = 'low'

                ot_row = (f'{qname}\t{n_genes}\t[{gids_s}]\t[{gnames_s}]\t[{cigar_s}]'
                          f'\t[{tnames_s}]\t[{ttypes_s}]'
                          f'\t{probe_gene}\t[{ot_gene_names_s}]\t[{ot_gene_types_s}]\t{concern}\n')
                fh_ot.write(ot_row)
    return no_hit

def main(args) -> None:
    print(message(f"TRACK module is aligning flipped probes to source transcripts", Mtype.START))
    # print(message(f"aligning query probes to target transcripts", Mtype.INFO))
    
    # if in all mode, make query fwd_oriented.fa
    if args.module == 'all':
        print(message(f"using fwd_oriented.fa probes", Mtype.INFO))
        args.query = os.path.join(args.out_dir, 'fwd_oriented.fa')

    if args.one_mismatch:
        afn = align_nm(args.query, args.target, "track", args)
    else:
        afn = align(args.query, args.target, "track", True, args)
    qfa = pyfastx.Fasta(args.query)

    # print(message(f"loading target transcriptome infos", Mtype.INFO))

    # if all mode, read in flip_t2g.csv else make a t2g
    if args.module == 'all':
        tinfos = load_tinfos(os.path.join(args.out_dir, 'flip_t2g.csv'))
        print(message(f"load in transcript to gene mappings from FLIP module", Mtype.INFO))
    else:
        fn = os.path.join(args.out_dir, 'track_t2g.csv')
        if not os.path.exists(fn) or args.force:
            # att_sep = ' ' if args.gtf else '='
            att_sep = check_annotation_ext(args.annotation)
            print(message(f"building transcript to gene mappings", Mtype.INFO))
            tinfos = build_tinfos(args.annotation, att_sep, args.schema, args.keep_dot)
            write_tinfos(fn, tinfos)
        else:
            tinfos = load_tinfos(fn)

    # print(message(f"detecting potential off-target probe activities", Mtype.INFO))
    if not args.one_mismatch:
        ainfos = track_target_pad(afn, qfa, args.pad_length, tinfos, not args.bowtie2,
                                  getattr(args, 'max_mismatches', -1))
        unaligned = get_unaligned(qfa, ainfos)
    else:
        tfa = pyfastx.Fasta(args.target)
        ainfos = track_target_one_mismatch(afn, qfa, tfa, tinfos)
        unaligned = get_unaligned(qfa, ainfos)
    print(message(f"{len(unaligned)} / {len(qfa)} probes unmapped. See track.unmapped.txt", Mtype.RESULT))
    write_lst2file(unaligned, os.path.join(args.out_dir, 'track.unmapped.txt'))
    no_hit = write_results(ainfos, args.out_dir)
    write_lst2file(no_hit, os.path.join(args.out_dir, 'track.no_hit.txt'))
    print(message(f"{len(no_hit)} / {len(qfa) - len(unaligned)} mapped probes with no passing hit. See track.no_hit.txt", Mtype.RESULT))
    print(message(f"finished tracking probes", Mtype.DONE))
