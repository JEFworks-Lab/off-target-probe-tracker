

# import libraries
import pandas as pd
from collections import defaultdict, Counter


# function to compute the difference: (total counts) - (self-match count)
def compute_diff(row):
    matches = row['Matches']
    total = sum(matches.values())
    # Get the count for the self-match if present; otherwise use 0
    self_count = matches.get(row['Gene'], 0)
    return total - self_count


# final df
bps_df = pd.DataFrame()

# iterate over the 16 bp files
for i in range(0, 16):

    # read the gene_summary file
    gene_summary = pd.read_csv(f"/home/caleb/Desktop/off-target-probe-checker/caleb/bowtie2/probe2targets.named.{0}bp.tsv", sep="\t", header=None)
    gene_summary.columns = ["probe_id", "n_targets", "target_names"]

    # get gene names
    gene_summary["probe_gene"] = gene_summary["probe_id"].str.split("|").str[1]
    # change column names
    gene_summary = gene_summary.rename(columns={"target_names": "probe_hit"})

    # create a defaultdict that will hold counter for each probe_gene
    gene_dict = defaultdict(Counter)

    # iterate over each row in the gene_summary dataframe
    for idx, row in gene_summary.iterrows():
        # key is target gene
        key = row['probe_gene']
        
        # get the probe_hit string and remove square brackets
        hits_str = row['probe_hit'].strip()  # remove leading/trailing whitespace
        if hits_str.startswith('[') and hits_str.endswith(']'):
            hits_str = hits_str[1:-1]
        
        # split the string on commas and strip any whitespace.
        hits = [hit.strip() for hit in hits_str.split(',') if hit.strip()]
        
        # add to counter for this probe_gene with each hit
        for hit in hits:
            gene_dict[key][hit] += 1

    # convert the defaultdict to a normal dict
    gene_dict = dict(gene_dict)

    gene_summary = pd.DataFrame(gene_dict.items(), columns=["Gene", "Matches"])

    # compute difference
    gene_summary['diff_' + str(i)] = gene_summary.apply(compute_diff, axis=1)

    # gene_summary.columns = ['Gene', 'Matches_' + str(i), 'diff_' + str(i)]/

    # merge with bps_df
    # if bps_df.empty:
    #     bps_df = gene_summary[['Gene', 'diff_' + str(i), 'Matches_' + str(i)]]
    # else:
    #     bps_df = bps_df.merge(gene_summary[['Gene', 'diff_' + str(i), 'Matches_' + str(i)]], on='Gene', how='outer')
    if bps_df.empty:
        bps_df = gene_summary[['Gene', 'diff_' + str(i)]]
    else:
        bps_df = bps_df.merge(gene_summary[['Gene', 'diff_' + str(i)]], on='Gene', how='outer')

print(bps_df)

# get rid of all rows that have 0 in all columns except the first column, that means no off-targets were detected ever
bps_df_filtered = bps_df.loc[~(bps_df.iloc[:, 1:] == 0).all(axis=1)]

# some genes have a different naming scheme, like NARS should be NARS1, 
# so remove all rows that have the same number all throughout their columns
# NOTE: this may get rid of ones i dont mean to, need to check naming schemes manually or change them
bps_df_filtered = bps_df_filtered.loc[~(bps_df_filtered.iloc[:, 1:].nunique(axis=1) == 1)]

bps_df_filtered


# plot the data
import matplotlib.pyplot as plt

# plot the where column diff is the x axis and the y axis is a line for each gene
for idx, row in bps_df_filtered.iterrows():
    plt.plot(row[1:], label=row['Gene'], marker='o')
    plt.xlabel('total bp difference')
    plt.ylabel('number of off-targets')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # rotate x axis labels
    plt.xticks(rotation=45)


