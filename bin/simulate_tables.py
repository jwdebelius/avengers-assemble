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

from qiime2 import Artifact


def define_regional_table(samples, seqs, read_length):
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
        lambda x: len(regex.findall('[RYSWKMBDHVN]', x))
        )
    len_check = seqs.astype(str).apply(lambda x: (len(x) >= read_length))
    degen_keep = degen_count.index[(degen_count <= 10) & len_check].values
    var_ids = [id_ for id_ in samples.index if (id_ in degen_keep)]
    samples = samples.loc[var_ids]
    table_seqs = seqs.copy().loc[var_ids]

    def _expand_degen(id_, seq_):
        degen_count = np.sum(seq_.degenerates())
        if degen_count == 0:
            return pd.Series({id_: seq_})
        else:
            return pd.Series({
                '%s@%s' % (id_, str(i).zfill(3)): exp 
                for i, exp in enumerate(seq_.expand_degenerates())
            })
    expand_seqs = pd.concat([_expand_degen(id_, seq) 
                             for id_, seq in table_seqs.iteritems()])

    # Checks the degeneate counts. Sequences with more than 10 degenerates 
    # get ignored because I have no clue how ot handle them?
    degen_ids = pd.Series({id_: id_.split('@')[0]
                           for id_ in expand_seqs.index}).reset_index()
    exp_ = degen_ids.groupby(0)['index'].apply(lambda x: x.values)
    exp_ = exp_.to_dict()

    # Does the sample subsampling and amplifies the table. We need to handle the 
    # dgenerate sequences seperately so that they're randomly selected and
    # subsampled
    samples = samples.astype(int)
    amp_table = []
    for sample in samples.columns:
        expanded = np.hstack([
            np.random.choice(exp_[id_], size=count, replace=True)
            for id_, count in 
            samples.loc[exp_.keys(), sample].iteritems()
            ])
        sub = pd.Series(np.random.choice(expanded, size=10000, replace=False))
        amp_table.append(sub.value_counts())

    amp_table = pd.concat(axis=1, objs=amp_table, sort=False).fillna(0)
    amp_table.columns = samples.columns

    exp_seqs = expand_seqs.loc[amp_table.index]

    return amp_table, exp_seqs.astype(str)


def hash_seq(x):
    """
    A lazy function to hash things
    """
    return hashlib.md5(x.encode()).hexdigest()


def define_asv_table(sample_table, sample_seqs, seq_length=260):
    """
    Makes an ASV table out of a regional table
    """
    # We get the subset of sequences and then group the table by these
    sample_seqs = \
        sample_seqs.loc[sample_seqs.apply(lambda x: len(x) >= seq_length)]
    asv_seqs = sample_seqs.apply(lambda x: x[:seq_length])
    sample_table['sequence'] = asv_seqs
    asv_table = sample_table.groupby('sequence').sum()    
    
    asv_biom = biom.Table(
        data=asv_table.values, 
        sample_ids=asv_table.columns,
        observation_ids=asv_table.index,
        )
    asv_seqs.index = asv_seqs.values
    asv_seqs.drop_duplicates(inplace=True)

    if set(asv_seqs) != set(asv_biom.ids(axis='observation')):
        raise ValueError("Sequence and table mismatch!")

    return asv_biom, asv_seqs


def _pcr_error(sample, sample_table, sample_seqs, rounds, pcr_error, read_length):
    """
    Simulates PCR error
    """

    # Sets up the sequence counts
    counts = np.hstack([
        [id_] * int(count) for id_, count in 
        sample_table[sample].dropna().iteritems()
        ])
    # Pulls out the sequences and converts them to integers
    sub_seqs = sample_seqs.loc[counts].copy().apply(lambda x: pd.Series(list(x)))
    sub_seqs = sub_seqs[np.arange(0, read_length)]
    sub_seqs.dropna(inplace=True)
    sub_seqs.replace({"A": 0, "C": 1, "G": 2, 'T': 3}, inplace=True)
    # Gets the error and random substitution
    error = (np.random.binomial(n=1, p=1e-4, size=sub_seqs.shape) * 
             np.random.randint(1, 3, size=sub_seqs.shape))
    # Adds the error and goes back to the sequence
    new_seq = ((sub_seqs.astype(int) + error) % 4)
    new_seq.replace({0: "A", 1: 'C', 2: "G", 3: "T"}, inplace=True)

    # Produces the dereplicated table
    return new_seq.apply(lambda x: ''.join(x), axis=1).value_counts()



def build_dereplicated_table(sample_table, sample_seqs, seq_length=260, 
                             error=1e-3, rounds=1):
    """
    Builds a deprelicated table
    """
    
    # Builds the dereplicated table and filters down to matching sequences
    sample_seqs = sample_seqs[np.arange(0, seq_length)].dropna()
    seq_ids = list(set(sample_seqs.index) & set(sample_table.index))
    sample_table = sample_table.loc[seq_ids].copy()
    sample_seqs = sample_seqs.loc[seq_ids]
    
    # Builds the dereplicated table
    derep_table = pd.concat(axis=1, sort=False, objs=[
        _pcr_error(sample, sample_table, sample_seqs, rounds, error, seq_length) 
        for sample in sample_table.columns
    ]).fillna(0)
    
    # Converts to a biom table
    derep_biom = biom.Table(
        data=derep_table.values,
        observation_ids=derep_table.index,
        sample_ids=sample_table.columns,
    )
    return derep_biom, pd.Series({seq_: seq_ for seq_ in derep_table.index})


def wrap_everything(seqs, table, read_length):
    """
    A useful wrapper function to keep everything together (hopefully)
    """
        
    # Builds the regional table
    region_table, region_seqs = define_regional_table(table, seqs, read_length)
    region_table.fillna(0, inplace=True)
    
    # Builds the ASV table
    asv_biom, asv_seqs = define_asv_table(region_table.copy(), 
                                          region_seqs.copy(), 
                                          read_length)
    asv_table = Artifact.import_data('FeatureTable[Frequency]', 
                                      asv_biom, biom.Table)
    asv_rep_seq = Artifact.import_data('FeatureData[Sequence]', 
                                       asv_seqs.dropna().apply(lambda x: skbio.DNA(x)))
    
    
    # Builds the dereplicated OTU set
    # Builds the dereplicated biom table
    derep_biom, derep_seqs = build_dereplicated_table(region_table.copy(), 
                                                      region_seqs.copy(), 
                                                      read_length)
    # Performs OTU clustering at 97% against gg 13-8
    derep_table = Artifact.import_data('FeatureTable[Frequency]', 
                                       derep_biom, biom.Table)

    derep_seqs = Artifact.import_data('FeatureData[Sequence]', 
                                      derep_seqs.apply(lambda x: skbio.DNA(x)))

    return asv_table, asv_rep_seq, derep_table, derep_seqs


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
    '--asv-table',
    help=('Path to save the ASV table')
    )
parser.add_argument(
    '--asv-rep-set',
    help=('Path to save the ASV representative sequences')
    )
parser.add_argument(
    '--dereplicated-table',
    help=("Location for the dereplicated table for OTU clustering")
    )
parser.add_argument(
    '--dereplicated-rep-set',
    help=("Location for the dereplicated sequences for OTU clustering")
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
    asv_table, asv_seqs, derep_table, derep_seqs = \
        wrap_everything(seqs, table, args.read_length)

    asv_table.save(args.asv_table)
    asv_seqs.save(args.asv_rep_set)
    derep_table.save(args.dereplicated_table)
    derep_seqs.save(args.dereplicated_rep_set)




