#!/usr/bin/env python

description = """
This script will build a set of synthetic samples based on resampling of the
american gut data. We'll do this by reading in the biom table and then 
randomly selecting a subset of samples using non ambigious samples. 

The samples will then be split into two approximately equal sized groups and 
used to generate a set of normalized, simulated samples
"""

import argparse
import os

import biom
import numpy as np
import pandas as pd

from qiime2 import Artifact


def subsample_data(metadata, table, age_split_fun, id_prefixes,
                   colors=None,
                   group_size=10, 
                   n_samples=20, seed=1234):
    """
    Generates a set of synthetic sequences based on the American Gut using
    seperation between people in their 20s and their 60s

    Parameters
    ----------
    keep_map : DataFrame
        A map of the samples with the age encoded in `age_cat` and the 
        index as the sample number
    samples : biom.Table
        The samples used for simulation
    age_split_fun: function
    group_size : int
        The number of samples to be combined to make the simulated sample
    split_size : int
        The number of samples to use in each set. (This is used to adjust)
        for uneven group sizes in the two age groups
    n_samples : int
        The number of samples to simulate for each age group

    Return 
    ------
    biom.Table
        All of the synthesized samples

    """
    # Sets a seed for consistency
    np.random.seed(seed)

    # Gets the age catgories 
    metadata['age_cat'] = metadata['age'].apply(age_split_fun)
    metadata.dropna(subset=['age_cat'], inplace=True)

    table.filter(metadata.index, axis='sample', inplace=True)
    table.add_metadata(metadata[['age_cat']].to_dict(orient='index'), axis='sample')

    print(id_prefixes)

    tables = [
        table.copy().filter(lambda v, id_, md: md['age_cat'] == prefix,
                            axis='sample', inplace=False)
        for prefix in id_prefixes
    ]

    synthetic = [
        tidy_names(pd.concat(axis=1, objs=[
            build_subsample(table, group_size, prefix) 
            for i in range(n_samples)]), 
        prefix)
        for prefix, table in zip(*(id_prefixes, tables))
    ]
    print(synthetic[0].head())

    all_samples = pd.concat(axis=1, sort=False, objs=synthetic)
    all_samples.fillna(0, inplace=True)

    all_biom = biom.Table(data=all_samples.values, 
                          sample_ids=all_samples.columns,
                          observation_ids=all_samples.index
                          )
    metadata = pd.DataFrame.from_dict(orient='index', data={
        id_: {'set': id_.split('.')[1],
              'age': id_.split('.')[2],
              }
        for id_ in all_biom.ids(axis='sample')
        })
    metadata.index.set_names('sample-id', inplace=True)
    metadata['composite'] = metadata.apply(lambda x: "{age}-{set}".format(**x),
                                            axis=1)
    if colors is not None:
        metadata['color'] = metadata['composite'].replace(colors)

    return all_biom, metadata

def tidy_names(table, prefix):
    table.rename(inplace=True, columns={
        int(i): f"sample.1.{prefix}.{i}" for i in table.columns
        })
    return table

def build_subsample(table, sample_size, group):
    """
    Builds a subsampled table from a biom object
    """
    sub_t = table.copy().subsample(sample_size, axis='sample', by_id=True)
    sub_t.filter(lambda v, id_, md: (v > 0).sum() > np.round(0.4 * sample_size), 
                 axis='observation', 
                 inplace=True)
    single = sub_t.copy().collapse(lambda id_, md: md['age_cat'], 
                                   axis='sample')
    single.norm(axis='sample', inplace=True)
    single.filter(lambda v, id_, md: (v > 2e-5).all(), 
                  axis='observation', 
                  inplace=True)
    single.norm(axis='sample', inplace=True)
    return pd.Series(np.round(single.data(group) * 5e5).round(0),
                     index=single.ids(axis='observation'))


parser = argparse.ArgumentParser(description=description)
parser.add_argument(
    '--keep-map',
    help=("The map of the samples to be used for simulation. This must "
          "contain  a sample identifier labeled `sample_name` and the age"
          "labeled age_cat"
          ),
    )
parser.add_argument(
    '--sample-table',
    help=("A biom table containing sample counts. The sample list should "
          "match those present in `--keep-map`"),
    )
parser.add_argument(
    '--simulated-table',
    help=('Filepath to the qiime artifact containing all the samples'),
    )
parser.add_argument(
    '--simulated-metadata',
    help=('Filepath to the qiime artifact containing all the samples'),
    )
parser.add_argument(
    '--samples-into-simulation',
    help=("The number of samples to combine into each simulated sample"),
    default=10,
    type=int,
    )
parser.add_argument(
    '--group-size',
    help=("The number of samples in each group to simulate for the set."),
    default=20,
    type=int,
    )
parser.add_argument(
    '--seed',
    help=('a random seed'),
    default=1234,
    type=int,
    )

if __name__ == '__main__':
    args = parser.parse_args()

    # Loads group
    keep_map = pd.read_csv(args.keep_map, sep='\t', dtype=str)   
    keep_map.set_index('sample_name', inplace=True)
    keep_map['age'] = keep_map['age'].astype(float)

    # Loads table
    sample_table = Artifact.load(args.sample_table).view(biom.Table)

    def cat_age(x):
        if pd.isnull(x):
            return x
        if x <= 3:
            return 'infant'
        elif (3 < x) & (x <= 12):
            return 'child'
        elif (20 < x) & (x <= 60):
            return 'adult'
        else:
            return np.nan

    colors = {'infant-1': '#1f78b4', 'infant-2': '#a6cee3', 
              'adult-1': '#e31a1c', 'adult-2': '#fb9a99'}


    all_biom, metadata = \
        subsample_data(metadata=keep_map.copy(), 
                       table=sample_table, 
                       age_split_fun=cat_age, 
                       id_prefixes=['infant', 'adult'],
                       group_size=args.samples_into_simulation,
                       n_samples=args.group_size,
                       seed=args.seed,
                       )
    all_table = Artifact.import_data('FeatureTable[Frequency]', all_biom)
    all_table.save(args.simulated_table)

    metadata.to_csv(args.simulated_metadata, sep='\t')



