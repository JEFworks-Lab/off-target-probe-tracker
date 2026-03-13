#!/usr/bin/env python

from opt.commons import *
from opt import flip, track, stat

VERSION="v0.0.3"

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--version', action='version', version=f'%(prog)s {VERSION}')

    # common args
    parser.add_argument('-o', '--out-dir', type=str, required=True, help="output dir")
    parser.add_argument('-p', '--threads', type=int, default=1, required=False, help="number of threads (default: 1)")
    parser.add_argument('--bam', action='store_true', default=False, required=False, help="")
    parser.add_argument('-b', '--binary', type=str, default=None, required=False, help="")
    parser.add_argument('--bowtie2', action='store_true', default=False, required=False, help="")
    # parser.add_argument('--gtf', action='store_true', default=False, required=False, help="") # deprecated
    parser.add_argument('-l', '--min-exact-match', required=False, help="", \
                default=20, type=int)
    parser.add_argument('--schema', type=lambda s: s.split(','), required=False, \
                default=['transcript', 'ID', 'Parent', 'gene_name', 'transcript_type'])
    parser.add_argument('--keep-dot', required=False, default=False, help="", \
                action='store_true')
    parser.add_argument('--force', required=False, help="", \
                default=False, action='store_true')
    parser.add_argument('--skip-index', required=False, help="", \
                default=False, action='store_true')
    
    subparsers = parser.add_subparsers(dest='module', \
                            help="[flip, track, stat, all]")

    # flip module
    parser_flip = subparsers.add_parser('flip', help="")
    parser_flip.add_argument('-q', '--query', type=str, required=True, \
                            help="query probe sequences (fasta)")
    parser_flip.add_argument('-t', '--target', type=str, required=True, \
                            help="target transcript sequences (fasta)")
    parser_flip.add_argument('-a', '--annotation', type=str, required=True, \
                    help="target transcriptome annotation (gff or gtf)")
    
    # track module
    parser_track = subparsers.add_parser('track', help="")
    parser_track.add_argument('-q', '--query', type=str, required=True, \
                            help="query probe sequences (fasta)")
    parser_track.add_argument('-t', '--target', type=str, required=True, \
                            help="target transcript sequences (fasta)")
    parser_track.add_argument('-a', '--annotation', type=str, required=True, \
                    help="target transcriptome annotation (gff or gtf)")
    parser_track.add_argument('-pl', '--pad-length', type=int, required=False, \
                    help="", default=0)
    parser_track.add_argument('-1', '--one-mismatch', action='store_true', \
                    default=False, required=False, help="")
    parser_track.add_argument('-mm', '--max-mismatches', type=int, default=-1, required=False,
                    help="Max mismatches allowed anywhere in the full probe (-1 = disabled)")
    
    # stat module
    parser_stat = subparsers.add_parser('stat', help="")
    parser_stat.add_argument('-i', '--in-file', type=str, required=True, \
                    help="track module results (i.e., probe2targets.csv)")
    parser_stat.add_argument('-q', '--query', type=str, required=True, \
                    help="query probe sequences (fasta)")
    parser_stat.add_argument('--exclude-pseudo', required=False, default=False, help="", \
                    action='store_true')
    parser_stat.add_argument('--pc-only', required=False, default=False, help="", \
                    action='store_true')
    parser_stat.add_argument('-s', '--syn-file', type=str, required=False,
                    help="", default=None)
    
    # all module
    parser_all = subparsers.add_parser('all', help="")
    parser_all.add_argument('-q', '--query', type=str, required=True, \
                            help="query probe sequences (fasta)")
    parser_all.add_argument('-t', '--target', type=str, required=True, \
                            help="target transcript sequences (fasta)")
    parser_all.add_argument('-a', '--annotation', type=str, required=True, \
                    help="target transcriptome annotation (gff or gtf)")
    parser_all.add_argument('-pl', '--pad-length', type=int, required=False, \
                            help="", default=0)
    parser_all.add_argument('-1', '--one-mismatch', action='store_true', \
                            default=False, required=False, help="")
    parser_all.add_argument('-mm', '--max-mismatches', type=int, default=-1, required=False,
                    help="Max mismatches allowed anywhere in the full probe (-1 = disabled)")
    parser_all.add_argument('--exclude-pseudo', required=False, default=False, help="", \
                    action='store_true')
    parser_all.add_argument('--pc-only', required=False, default=False, help="", \
                    action='store_true')
    parser_all.add_argument('-s', '--syn-file', type=str, required=False,
                    help="gene synonyms CSV file with 2 columns", default=None)
    parser_all.add_argument('-i', '--in_file',
                    help="path to input file (default: probe2targets.tsv in out_dir)",
                    default=None)
    
    args = parser.parse_args()
    if args.module not in ['flip', 'track', 'stat', 'all']:
        parser.error(f"Invalid module {args.module}. Valid options are: flip, track, predict, all")

    # For the 'all' module, handle gene synonyms if provided
    if args.module == 'all':
        if args.syn_file is not None:

            # directory where the original query FASTA lives
            query_dir = os.path.dirname(args.query)
            if query_dir == "":
                query_dir = "."

            # filename & extension
            base_name = os.path.basename(args.query)
            base, ext = os.path.splitext(base_name)
            if ext == "":
                ext = ".fa"

            # save the updated FASTA
            args.query_genesyn = os.path.join(query_dir, f"{base}_genesyn{ext}")

            print(message(f"writing synonym-updated query FASTA to {args.query_genesyn}", Mtype.INFO))
            replace_genenames_with_synonyms(
                query_fasta=args.query,
                gene_syn_csv=args.syn_file,
                out_fasta=args.query_genesyn
            )

        else:
            args.query_genesyn = args.query

    return args


# confirm directory exists, if not create
def check_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

# check if flip args are present
def check_flip_args(args) -> bool:
    check_dir(args.out_dir)
    return all(os.path.exists(pth) for pth in \
            [args.query, args.target, args.annotation])

# check if track args are present
def check_track_args(args) -> bool:
    check_dir(args.out_dir)
    return all(os.path.exists(pth) for pth in \
            [args.query, args.target, args.annotation])

# check if stat args are present
def check_stat_args(args) -> bool:
    check_dir(args.out_dir)
    return all(os.path.exists(pth) for pth in \
            [args.in_file])

# check if all args are present
def check_all_args(args) -> bool:
    check_dir(args.out_dir)
    return all(os.path.exists(pth) for pth in \
            [args.query, args.target, args.annotation])

def main() -> None:
    args = parse()

    if args.module == 'flip':
        if not check_flip_args(args):
            print(message(f"cannot locate files", Mtype.ERROR))
            sys.exit(-1)
        print(message(f"### FLIP ###", Mtype.START))
        param_fn = os.path.join(args.out_dir, "flip_params.json")
        store_params(args, param_fn)
        flip.main(args)
    elif args.module == 'track':
        if not check_track_args(args):
            print(message(f"cannot locate files", Mtype.ERROR))
            sys.exit(-1)
        if args.one_mismatch and args.pad_length > 0:
            print(message(f"cannot use -1 mode with positive pad_length", Mtype.ERROR))
            sys.exit(-1)
        if args.one_mismatch and args.max_mismatches >= 0:
            print(message(f"cannot use --one-mismatch with --max-mismatches", Mtype.ERROR))
            sys.exit(-1)
        print(message(f"### TRACK ###", Mtype.START))
        param_fn = os.path.join(args.out_dir, "track_params.json")
        store_params(args, param_fn)
        track.main(args)
    elif args.module == 'stat':
        if not check_stat_args(args):
            print(message(f"cannot locate files", Mtype.ERROR))
            sys.exit(-1)
        if args.pc_only and args.exclude_pseudo:
            print(message(f"cannot use both --pc-only and --exclude-pseudo flags", Mtype.ERROR))
            sys.exit(-1)
        print(message(f"### STAT ###", Mtype.START))
        param_fn = os.path.join(args.out_dir, "stat_params.json")
        store_params(args, param_fn)
        stat.main(args)
    # adding an option to do all three modules in one go
    elif args.module == 'all':
        # if a synonym FASTA was created, use it as the query for downstream modules
        if hasattr(args, "query_genesyn") and args.query_genesyn is not None:
            args.query = args.query_genesyn
        # confirm we have files and nothing overlaps
        if not check_all_args(args):
            print(message(f"cannot locate files", Mtype.ERROR))
            sys.exit(-1)
        if args.one_mismatch and args.pad_length > 0:
            print(message(f"cannot use -1 mode with positive pad_length", Mtype.ERROR))
            sys.exit(-1)
        if args.one_mismatch and args.max_mismatches >= 0:
            print(message(f"cannot use --one-mismatch with --max-mismatches", Mtype.ERROR))
            sys.exit(-1)
        if args.pc_only and args.exclude_pseudo:
            print(message(f"cannot use both --pc-only and --exclude-pseudo flags", Mtype.ERROR))
            sys.exit(-1)
        print(message(f"### running the ALL module ###", Mtype.START))

        # print out param json
        param_fn = os.path.join(args.out_dir, "all_params.json")
        # store params
        store_params(args, param_fn)
        # run opt modules
        flip.main(args)
        track.main(args)
        stat.main(args)

if __name__ == "__main__":
    main()
