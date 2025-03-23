from otpc.commons import *

def convert_md2bit(s):
    running = ""
    bit_s = ""
    for c in s:
        if c.isdigit():
            running += c
        else:
            if len(running) > 0:
                bit_s += '1' * int(running)
            bit_s += '0'
            running = ""
    if len(running) > 0:
        bit_s += '1' * int(running)
    return bit_s

def convert_md2bit_nucmer(s, tstart):
    running = ""
    bit_s = ""
    for c in s:
        if c.isdigit():
            running += c
        else:
            if len(running) > 0:
                bit_s += '1' * int(running)
            bit_s += '0'
            running = ""
    if len(running) > 0:
        bit_s += '1' * (int(running) - tstart)
    return bit_s

def convert_md2bit_del(s):
    running = ""
    bit_s = ""
    ignore = False
    mismatch_info = []
    for c in s:
        if c.isdigit():
            ignore = False
            running += c
        else:
            if ignore: continue
            if c == '^':
                if len(running) > 0:
                    bit_s += '1' * int(running)
                running = ""
                ignore = True
                continue
            else:
                if len(running) > 0:
                    bit_s += '1' * int(running)
                mismatch_info.append(len(bit_s))
                bit_s += '0'
                running = ""
    if len(running) > 0:
        bit_s += '1' * int(running)
    return bit_s, mismatch_info

def convert_md2bit_nucmer_del(s, tstart):
    running = ""
    bit_s = ""
    ignore = False
    mismatch_info = []
    for c in s:
        if c.isdigit():
            ignore = False
            running += c
        else:
            if ignore: continue
            if c == '^':
                if len(running) > 0:
                    bit_s += '1' * int(running)
                running = ""
                ignore = True
                continue
            else:
                if len(running) > 0:
                    bit_s += '1' * int(running)
                mismatch_info.append(len(bit_s))
                bit_s += '0'
                running = ""
    if len(running) > 0:
        bit_s += '1' * (int(running) - tstart)
    return bit_s, mismatch_info

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
    out = ""
    assert len(s1) == len(s2)
    for i in range(len(s1)):
        if s1[i] == '1' and s2[i] == '1':
            out += '1'
        else:
            out += '0'
    return out

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

def track_target_pad(fn, qfa, pad, tinfos, is_nucmer) -> dict:
    ainfos = dict()
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
                qlen = len(qfa[qname].seq)
                crit_bvec = "0" * pad + "1" * (qlen - 2 * pad) + "0" * pad
                assert len(crit_bvec) == qlen # sanity check
                crit_dvec = int(crit_bvec, 2)
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
                        if crit_dvec & int(md_bvec, 2) == crit_dvec:
                            ainfos[qname].add((tname, (tinfos[tname][0], tinfos[tname][1]), \
                                        tinfos[tname][2], compress_bvec(final_bvec)))
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
                        if hit:
                            assert final_bvec is not None
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
                        if crit_dvec & int(final_bvec, 2) == crit_dvec:
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

def check_lft_and_rgt(mrec, qname, qry_fa, tgt_fa, max_nm):
    qlen = len(qry_fa[qname].seq)
    tname, tst, qst, mlen = mrec
    if mlen == 40:
        return (True, '1' * mlen, 0)
    qst -= 1
    tst -= 1
    qen = qst + mlen
    ten = tst + mlen
    lft_qos = qst
    rgt_qos = qlen - qen
    qseq = qry_fa[qname].seq
    tseq = tgt_fa[tname].seq
    lft_tseq = tseq[tst - lft_qos:tst]
    rgt_tseq = tseq[ten:ten + rgt_qos]
    lft_qseq = qseq[qst - lft_qos:qst]
    rgt_qseq = qseq[qen:qen + rgt_qos]
    if len(lft_tseq) != lft_qos: # tseq runs out at 5'
        return (False, None, -1)
    if len(rgt_tseq) != rgt_qos: # tseq runs out at 3'
        return (False, None, -1)
    lft_mvec = ""
    lft_nm = 0
    for i in range(lft_qos):
        if lft_tseq[i] == lft_qseq[i]:
            lft_mvec += "1"
        else:
            lft_mvec += "0"
            lft_nm += 1
    rgt_mvec = ""
    rgt_nm = 0
    for i in range(rgt_qos):
        if rgt_tseq[i] == rgt_qseq[i]:
            rgt_mvec += "1"
        else:
            rgt_mvec += "0"
            rgt_nm += 1
    mvec = lft_mvec + ('1' * mlen) + rgt_mvec
    assert len(mvec) ==  qlen # sanity check
    nm = lft_nm + rgt_nm
    return (nm <= max_nm, mvec, nm)

def track_target_nm(fn, qfa, tfa, max_nm, tinfos) -> dict:
    mums = load_mums(fn)
    ainfos = dict()
    for qname in mums:
        ainfos[qname] = set()
        for mrec in mums[qname]:
            tname = mrec[0]
            is_pass, mvec, _ = check_lft_and_rgt(mrec, qname, qfa, tfa, max_nm)
            if is_pass:
                ainfos[qname].add((tname, (tinfos[tname][0], tinfos[tname][1]), \
                                            tinfos[tname][2], compress_bvec(mvec)))
    return ainfos
    
def write_results(ainfos, d) -> list:
    fn = os.path.join(d, 'probe2targets.tsv')
    no_hit = []
    with open(fn, 'w') as fh:
        fh.write('probe_id\tn_genes\tgene_ids\tgene_names\tcigars\ttranscript_ids\ttranscript_types\n')
        for qname in ainfos:
            if len(ainfos[qname]) == 0:
                no_hit.append(qname)
                continue
            tnames = [x[0] for x in ainfos[qname]]
            genes = [x[1] for x in ainfos[qname]]
            gids = [x[0] for x in genes]
            gnames = [x[1] for x in genes]
            ttypes = [x[2] for x in ainfos[qname]]
            cigars = [x[3] for x in ainfos[qname]]
            try:
                assert len(gids) == len(gnames) # sanity check
            except:
                print(message(f">1 reference gene IDs might share the same gene name", Mtype.WARN))
                print(gids)
                print(gnames)
            gids_s = ','.join(gids)
            # handle no gene names
            gnames_s = ','.join('None' if x is None else x for x in gnames)
            cigar_s = ','.join(cigars)
            ttypes_s = ','.join(ttypes)
            tnames_s = ','.join(tnames)
            # n_genes = # of distinct gene_names
            fh.write(f'{qname}\t{len(set(gnames))}\t[{gids_s}]\t[{gnames_s}]\t[{cigar_s}]\t[{tnames_s}]\t[{ttypes_s}]\n')
    return no_hit

def main(args) -> None:
    print(message(f"aligning query probes to target transcripts", Mtype.PROG))
    if args.one_mismatch:
        afn = align_nm(args.query, args.target, "track", args)
    else:
        afn = align(args.query, args.target, "track", True, args)
    qfa = pyfastx.Fasta(args.query)

    print(message(f"loading target transcriptome infos", Mtype.PROG))
    fn = os.path.join(args.out_dir, 'track_t2g.csv')
    if not os.path.exists(fn) or args.force:
        att_sep = ' ' if args.gtf else '='
        print(message(f"building t2g mappings", Mtype.PROG))
        tinfos = build_tinfos(args.annotation, att_sep, args.schema, args.keep_dot)
        write_tinfos(fn, tinfos)
    else:
        tinfos = load_tinfos(fn)
    
    print(message(f"detecting potential off-target probe activities", Mtype.PROG))
    if not args.one_mismatch:
        ainfos = track_target_pad(afn, qfa, args.pad_length, tinfos, not args.bowtie2)
        unaligned = get_unaligned(qfa, ainfos)
    else:
        tfa = pyfastx.Fasta(args.target)
        ainfos = track_target_nm(afn, qfa, tfa, args.max_mismatch, tinfos)
        unaligned = get_unaligned(qfa, ainfos)
    print(message(f"{len(unaligned)} / {len(qfa)} probes unmapped", Mtype.PROG))
    write_lst2file(unaligned, os.path.join(args.out_dir, 'track.unmapped.txt'))
    no_hit = write_results(ainfos, args.out_dir)
    write_lst2file(no_hit, os.path.join(args.out_dir, 'track.no_hit.txt'))
    print(message(f"{len(no_hit)} / {len(qfa) - len(unaligned)} mapped probes with no passing hit", Mtype.PROG))
    print(message(f"finished", Mtype.PROG))