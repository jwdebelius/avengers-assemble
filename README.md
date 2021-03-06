# Avengers Assemble

This repository contains the analysis workflow for Debelius et al (doi: [10.1101/2021.03.23.436606](https://www.biorxiv.org/content/10.1101/2021.03.23.436606v1)). 


## Installation

The analysis was done based on a qiime2-2020.11 environment with [RESRIPt](https://github.com/bokulich-lab/RESCRIPt) and [q2-sidle](https://q2-sidle.readthedocs.io/en/latest/) installed. 

To build this environment, install [qiime2-202011](https://docs.qiime2.org/2020.11/install/) according to your system requirements. Then, follow the [q2-sidle](https://docs.qiime2.org/2020.11/install/) installation instructions. You will also need snakemake > 5.3.

```bash
conda install --channel bioconda snakemake=5.3
```

Finally, you will need to clone and download the repository in whatever way works best for you. 


## Input Files

### Simulation

To run the simulation, you will need to download the 100nt OTU table (ID 45113, CRC32: 7066cf83) and metadata from [Qiita Study 850](https://qiita.ucsd.edu/study/description/850). The OTU table should be imported into QIIME 2 using the command,

```bash
qiime tools import \
 --type 'FeatureTable[Frequency]' \
 --source-fomat BIOMV210Format \
 --input-path 45113_otu_table.biom \
 --output-path start-table.qza

```

And placed in the folder specified as the config file (`input_dir` under `simulation`). (By default, this folder is `ipynb/data/inputs/simulation`).

The metadata should also be downloaded from the qiita study (found under the sample information tab), and renamed to `start-metadata.txt` and moved to the simulation input folder `ipynb/data/inputs/simulation`.

### Real Data

Real data can be downloaded from ENA study [PRJEB37382](https://www.ebi.ac.uk/ena/browser/view/PRJEB37382), accessions ERR4704801-ERR4704848. The should be placed in the location specified in the config file. By default, they are expected in the `ipynb/data/inputs/real/seqs` directory. Metadata is inferred from the filenames. 

### Database Files

The analysis is currently based on the Silva 128 QIIME release, available [here](https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_release.tgz). The files should be imported into QIIME 2 and saved to the `ipynb/data/reference/silva-ori/` folder. The SEPP reference tree can be downloaded from the [qiime2 resources page](https://data.qiime2.org/2021.2/common/sepp-refs-silva-128.qza). The formatted Optivag and Greengenes databases are already in the repository.


## Running the workflow

You can be through easily through snakemake based on the current repository set up.

```bash
snakemake 
```

## Workflow steps

### Simulation

1. **Reference Simulation**. The greengenes database is used to filter the reference OTU table to exclude sequences with more than 5 degenerates, and then a reference dataset is simulated from the original data by averaging abundance across multiple individuals of the same age group.
	* `filter_degenerete_reference`
	* `filter_ref_table`
	* `simulate_reference_samples`
2. **Regional Simulation**. Regions of the 16s gene are extracted and combined with the full length reference table to build a regional "ASV" table and representative sequence set. These are the starting points for our three test methods.
	* `extract_regions`
	* `generate_asvs_and_derep_tables`
3. **Database Preparation**. The database is prepared by filtering to remove sequences with a large number of degenerate nucleotides from the reference and database sequences which do not have class-level annotations. These are used to generate per-region 
	* `filter_silva_reference`
	* `extract_sidle_regions`
	* `prepare_extracted_region`
4. **Reference Reconstruction**. The simulated reference table, the greengenes 13_8 taxonomy and the greengenes 13_8 tree are copied from their respective locations.
	* `copy_sim_reference`
5. **OTU Clustering**. The table is generated by merging all ASV representative sequences and the feature tables from all the individual regions and then clustering them closed reference at 99% identity against the Silva 128 reference database 
	* `concatenate_single_asv_counts`
	* `concatenate_single_rep_seqs`
	* `cluster_otus`
	* `get_otu_tree`
	* `get_otu_taxonomy`
6. **ASV denosing**. The table is generated by merging all ASV counts. Taxonomy comes from naive bayesian classification of all sequences; sequences which don't have phylum level annotation are excluded. The tree is built by fragment insertion of the merged sequences
	* `concatenate_single_asv_counts`
	* `concatenate_single_rep_seqs`
	* `classify_taxonomy`
	* `get_asv_rep_seq`
	* `get_asv_table`
	* `build_insertion_tree`
6. **Sidle**. For Sidle, the regional sequences are aligned to the regional kmers. The alignments, regional database maps, and feature counts are used to reconstruct the table. Taxonomy is reconstructed using the database map to identify unresolved sequences. To build the phylogenetic tree, we generate consensus fragments and then insert them into them in the Silva 128 backbone.
	* `align_sidle_regions`
	* `reconstruct_table_single`
	* `reconstruct_taxonomy`
	* `reconstruct_fragments_single`
	* `reconstructed_fragment_insertion`

7. **Alpha Diversity**. Within sample (alpha) diversity is calculated from multiple rarefaction. We look at observed features, Shannon, Pielou's evenness and faith's diversity. The results are summarized using the `Table1-CompareAlpha.ipynb` notebook. This generated Table 1.
	* `rarefy_iteration`
	* `alpha_iter`
	* `alpha_table`

8. **Beta Diversity**. The between sample (beta) diversity is also calculated from rarified tables, and then, we look at three feature-based metrics (unweighted UniFrac, weighted UniFrac, and Bray-Curtis) as well as a collapsed metric (genus level Bray-Curtis). The results are summarized in the `ipynb/Table2-CompareBeta.ipynb` notebook. This generates Table 2.
	* `rarefy_iteration`
	* `genus_braycurtis`
	* `beta_phylogenetic`
	* `braycurtis`
	* `adonis`
	* `beta_table`

9. **Taxonomic Comparison**. We compared the taxonomy at class level using the merged, unnormalized tables through the `ipynb/Figure1-Taxonomy.ipynb`. This generates Figure 2 and Table S3.
	* `compare_taxonomy`

### Real Data Analysis

1. **Data Preparation**. We import the sequences into QIIME 2, and then denoise them using dada2 before denoising the sequences.
	* `build_manifest`
	* `import_seqs`
	* `trim_data_v13`
	* `denoise_v13`
	* `trim_data_v34`
	* `denoise_v34`
	* `trim_posthoc`

2. **Database Preparation**. For benchmarking the greengenes database needs to be prepared (`filter_greengenes_reference`), otherwise, the Optivag database is assumed to be tidied. Then, we extract the regions and prepare them from sidle alignment.
	* `filter_greengenes_reference`
	* `extract_sidle_regions`
	* `prepare_extracted_region`

3. **Reference Reconstruction**. The V34 region alone was used as the reference for the simulation. We used the paired end table, classified the taxonomy using a naive bayesian classifier for the full length sequences, and then filtered the table to remove anything missing a phylum level or kingdom level designation.
	* `classify_regional_taxonomy`
	* `get_real_reference_data`
	* `get_real_representative_seqs`
5. **OTU Clustering**. The table is generated by merging all ASV representative sequences and the feature tables from all the individual regions and then clustering them closed reference at 99% identity
	* `concatenate_single_asv_counts`
	* `concatenate_single_rep_seqs`
	* `cluster_otus`
	* `get_otu_taxonomy`
6. **ASV denosing**. The table is generated by merging all ASV counts. Taxonomy comes from naive bayesian classification of all sequences; sequences which don't have phylum level annotation are excluded. The tree is built by fragment insertion of the merged sequences
	* `concatenate_single_asv_counts`
	* `concatenate_single_rep_seqs`
	* `classify_taxonomy`
	* `get_asv_rep_seq`
	* `get_asv_table`

6. **Sidle**. For Sidle, the regional sequences are aligned to the regional kmers. The alignments, regional database maps, and feature counts are used to reconstruct the table. Taxonomy is reconstructed using the database map to identify unresolved sequences.
	* `align_sidle_regions`
	* `reconstruct_table_single`
	* `reconstruct_taxonomy`

7. **Data Synthesis**. The taxonomic annotation was checked and the tables evaluated using the `ipynb/Figure2-RealData.ipynb` notebook.
	* `real_data_figure`

## Benchmarking

The snake file contains a section dedicated to preparing the sidle reads, and it will output a table. This will generate a series of benchmarking files, in the benchmark folder. The `Matlab` folder contains the modified m scripts to run with SMURF.


