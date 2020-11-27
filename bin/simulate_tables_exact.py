#!/usr/bin/env python

description = """
Simulates the an ASV table and dereplicated sequences for OTU clustering
and meta analysis exploration
"""


import argparse
import hashlib

import biom
import numpy as np
import pandas as pd
import regex
import skbio
from skbio import DNA

from qiime2 import Artifact


def build_regional_table(samples, seqs, read_length, depth=10000):
    """
    Amplifies a regional table from the original table

    Parameters
    ----------
    samples : DataFrame
        A biom-style dataframe where the index is sequences (corresponding)
        to the sequences in seqs

    """
    # Filters the sequences and table down to those present in the sequences
    degen_count = seqs.astype(str).apply(
        lambda x: len(regex.findall('[RYSWKMBDHVN]', str(x)))
        )
    len_check = seqs.astype(str).apply(lambda x: (len(x) >= np.absolute(read_length)))
    print(len_check.sum())
    degen_keep = degen_count.index[(degen_count <= 10) & len_check].values
    var_ids = [id_ for id_ in samples.index if (id_ in degen_keep)]
    samples = samples.loc[var_ids]

    table_seqs = seqs[var_ids].copy()

    # Expands the degenerate sequences
    exp_ = {id_: np.array([str(s) for s in seq_.expand_degenerates()])
            for id_, seq_ in table_seqs.iteritems()}
    
    # Generates the table by expanding the sequence where degenerates are
    # randomnly sampled and then combining the sequences into a table
    forward_reads = dict([])
    reverse_reads = dict([])
    joined_reads = dict([])
    for sample in samples.columns:
        expanded = pd.Series(np.hstack([
                np.random.choice(exp_[id_], size=count, replace=True)
                for id_, count in 
                samples.loc[samples[sample] > 0, sample].iteritems()])
        ).apply(lambda x: x[:read_length])
        joined_reads[sample] = expanded.value_counts()
 #         joined_reads[sample] = pd.Series(np.random.choice(expanded, depth, replace=False)
#                                          ).value_counts()
    joined_table = pd.DataFrame.from_dict(orient='index', data=joined_reads).fillna(0)
    joined_table = joined_table.loc[joined_table.sum(axis=1) > 10]
    joined_reads = pd.Series({
        hash_seq(seq): DNA(seq) 
        for seq in joined_table.columns
    })
    joined_table.rename(columns={seq: hash_seq(seq) for seq in joined_table.columns},
                        inplace=True)
    jointable = Artifact.import_data('FeatureTable[Frequency]', joined_table, pd.DataFrame)
    joinseqs = Artifact.import_data('FeatureData[Sequence]', joined_reads, pd.Series)
    
    return jointable, joinseqs


def hash_seq(x):
    """
    A lazy function to hash things
    """
    return hashlib.md5(x.encode()).hexdigest()




parser = argparse.ArgumentParser(description=description)
parser.add_argument(
    '--sequences',
    help=('The extracted regional sequences'),
    )
parser.add_argument(
    '--sample-table',
    help=("A table of the samples to be simulated."),
    )
parser.add_argument(
    '--read-length',
    help=("length of the final simulated sequences. We treat these as paired"
          "end sequences"),
    type=int,
    )
parser.add_argument(
    '--depth',
    help=('Depth of the final samples'),
    default=10000,
    )
parser.add_argument(
    '--table',
    help=('Path to save the joined-sequence table')
    )
parser.add_argument(
    '--reads',
    help=('Path to save the representative sequences for the joined table'),
    )
parser.add_argument(
    '--random',
    help=("A seed for randomn number generation"),
    default=None,
    )

if __name__ == '__main__':
    args = parser.parse_args()
    
    # Gets the regional sequences
    seqs =  Artifact.load(args.sequences).view(pd.Series)

    # Gets the table
    table = Artifact.load(args.sample_table).view(pd.DataFrame).fillna(0).astype(int).T

    # Simulates the samples
    if args.random is not None:
        np.random.seed(int(args.random))
    # fwd_table, rev_table, join_table, fwd_seqs, rev_seqs, join_seqs = \
    join_table, join_seqs = \
        build_regional_table(table, seqs, args.read_length, depth=args.depth)

    # fwd_table.save(args.fwd_table)
    # fwd_seqs.save(args.fwd_reads)
    # rev_table.save(args.rev_table)
    # rev_seqs.save(args.rev_reads)
    join_table.save(args.table)
    join_seqs.save(args.reads)




