

####################################################################################

# import libraries
import pandas as pd
from collections import Counter, defaultdict


####################################################################################


### read in the probe2targets.tsv file ###

# file_name
file_name = "out_6bp_nucmer_hubmap"

# read in the probe2targets.tsv file
probe2targets = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/otpc/{file_name}/probe2targets.tsv", sep="\t")


####################################################################################


### helper function ###

# parse bracketed lists
def parse_bracketed_list(value: str):
    if not value:
        return []
    # remove square brackets
    trimmed = value.strip("[]")
    if not trimmed.strip():
        return []
    # split on commas and strip whitespace/quotes from each item
    return [x.strip().strip("'\"") for x in trimmed.split(",") if x.strip().strip("'\"")]


####################################################################################


### creating df ###

# create a defaultdict that will hold counter for each probe_gene
gene_dict = defaultdict(
    lambda: defaultdict(
        lambda: {
            'count': 0,
            'cigars': list(),
            'transcript_types': list(),
            'target_names': list()
        }
    )
)


# go over the rows of the DataFrame
for idx, row in probe2targets.iterrows():
    # splot "probe_id" 
    parts = row['probe_id'].split('|')
    # if there are not 3 parts, skip NOTE: might need to change in future
    if len(parts) != 3:
        continue
    _, gene_name, probe_id = parts

    # Convert each bracketed field into a list of items
    cigars_list = parse_bracketed_list(str(row['cigars']))
    transcript_types_list = parse_bracketed_list(str(row['transcript_types']))
    target_names_list = parse_bracketed_list(str(row['target_names']))

    # Update dictionary
    gene_dict[gene_name][probe_id]['count'] = len(cigars_list)
    gene_dict[gene_name][probe_id]['cigars'] = cigars_list
    gene_dict[gene_name][probe_id]['transcript_types'] = transcript_types_list
    gene_dict[gene_name][probe_id]['target_names'] = target_names_list

# convert the internal sets to lists
final_dict = {}
for gene_name, probe_data in gene_dict.items():
    final_dict[gene_name] = {}
    for probe_id, details in probe_data.items():
        final_dict[gene_name][probe_id] = {
            'count': details['count'],
            'cigars': details['cigars'],
            'transcript_types': details['transcript_types'],
            'target_names': details['target_names'],
        }

# print(final_dict)


# aggregate the data
aggregator = defaultdict(dict)

for gene_name, probes in final_dict.items():
    for probe_id, details in probes.items():
        for tname in details['target_names']:
            # get index of the gene_name in the target_names list
            gene_idx = details['target_names'].index(tname)
            if tname not in aggregator[gene_name]:
                # if first time seeing this tname for this gene
                ls_tt = list()
                ls_cigars = list()
                ls_tt.append(details['transcript_types'][gene_idx])
                ls_cigars.append(details['cigars'][gene_idx])
                aggregator[gene_name][tname] = {
                    'count': 1,
                    'transcript_types': ls_tt,  # store as set or list
                    'cigars': ls_cigars
                }
            else:
                # weve seen this (gene_name, tname) before, so just update the count & cigars
                aggregator[gene_name][tname]['count'] += 1
                # aggregator[gene_name][tname]['cigars'].append(details['cigars'][gene_idx])
                # Do NOT update transcript_types here to avoid duplicates


# Now convert aggregator to a DataFrame
rows = []
for gene_name, target_dict in aggregator.items():
    aligned_to_list = []
    total_hits_list = []
    transcript_types_list = []
    cigars_list = []
    
    # Sort tnames for consistent ordering
    for tname in sorted(target_dict.keys()):
        info = target_dict[tname]
        aligned_to_list.append(tname)
        total_hits_list.append(str(info['count']))
        transcript_types_list.append(",".join(sorted(info['transcript_types'])))
        cigars_list.append(",".join(sorted(info['cigars'])))
    
    rows.append((
        gene_name,
        ",".join(aligned_to_list),
        ",".join(total_hits_list),
        ",".join(transcript_types_list),
        ",".join(cigars_list)
    ))

df_final = pd.DataFrame(rows, columns=["Gene", "aligned_to", "total_hits", "transcript_types", "cigars"])


# iterate over the rows of the DataFrame and remove rows that the only aligned_to gene is not the same as Gene
# aka getting genes with off-targets
gene_summary_filtered = df_final[
    df_final.apply(
        lambda row: len(row['aligned_to'].split(',')) > 1 or row['aligned_to'].split('(')[0] != row['Gene'],
        axis=1)]


# FINAL OUTPUT
gene_summary_filtered

