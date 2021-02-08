###############################################################################
#                            Gets reference dataset                           #
###############################################################################

import os
import numpy as np

configfile : 'config.yaml'
ruleorder:  collapse_to_level > rarify_table

rule all:
    input:
        "data/output/tables/table-1-alpha.tsv",
        "data/output/tables/table_2_beta_rarefaction.tsv",
        'data/output/figures/figure2/fig2_class_relative_abund.png',
        'data/output/tables/table_s2_class_comparison.tsv'
        'data/output/figures/figure3/figure2_vaginal.png'


###############################################################################
#                            Gets reference dataset                           #
###############################################################################
#                                                                             #
# Refrence data was downloaded from the global gut study in qiita (Study ID   #
# 805) and imported into QIIME 2.                                             #
#                                                                             #
# This block will take the reference sequences fitler them to anything with   #
# fewer than 5 degenerate nucloetides (`filter_degenerate_reference`,         #
# `filter_ref_table`). Then, it will simulate a set of samples based on the   #
# input data. We'll use infants (children under 3) and adults (ages 20 - 60)  #
# as the seeds for the analyis.                                               #
#                                                                             #
###############################################################################

rule filter_degenerate_reference:
  input:
    config['simulation']['taxonomy']['reference']['rep_seq']
  output:
    expand("{ref_dir}/gg_97_5degen.qza", ref_dir=config['simulation']['ref_dir'])
  shell:
    """
    qiime sidle filter-degenerate-sequences \
     --i-sequences {input} \
     --p-n-workers 2 \
     --p-max-degen 5 \
     --o-filtered-sequences {output}
    """

rule filter_ref_table:
  input:
    table=expand('{input_dir}/start-table.qza', input_dir=config['simulation']['input_dir']),
    sequences=expand("{ref_dir}/gg_97_5degen.qza", ref_dir=config['simulation']['ref_dir']),
  output:
    expand("{input_dir}/start-table_5degen.qza", input_dir=config['simulation']['input_dir'])

  shell:
    """
    qiime feature-table filter-features \
     --i-table {input.table} \
     --m-metadata-file {input.sequences} \
     --o-filtered-table {output}
    """

rule extract_reference_subset:
  input:
    table = expand("{input_dir}/start-table_5degen.qza", input_dir=config['simulation']['input_dir']),
    metadata = expand('{input_dir}/start-metadata.txt', input_dir=config['simulation']['input_dir']),
  output:
    table = '{outdir}/simulation/samples/table.qza', 
    metadata = '{outdir}/simulation/samples/metadata.tsv',
  params:
    group_size = config['simulation']['group_size']
  shell:
    """
    python bin/simulate_samples.py \
     --keep-map {input.metadata} \
     --sample-table {input.table} \
     --group-size {params.group_size} \
     --simulated-table {output.table} \
     --simulated-metadata {output.metadata}
    """


################################################################################
#                            Simulates Test Dataset                            #
################################################################################
#                                                                              #
# This block uses the simulated data to generate ASV tables based on the       # 
# reference samples in all the regions that will be used. We'll then cluster   #
# those sequences closed reference and perform a kmer-based alignment of the   #
# sequences for each dataset.                                                  #
#                                                                              #
################################################################################

rule extract_regions:
  input:
    lambda w: config['simulation']['taxonomy']['reference']['rep_seq']
  output:
    '{outdir}/simulation/asv_ref/{region}-seqs.qza'
  params:
    fwd = lambda w: config['simulation'][w.region]['fwd'],
    rev = lambda w: config['simulation'][w.region]['rev']
  shell:
    """
    qiime feature-classifier extract-reads \
     --i-sequences {input} \
     --p-f-primer {params.fwd} \
     --p-r-primer {params.rev} \
     --p-min-length 50 \
     --p-max-length 700 \
     --p-n-jobs 5 \
     --o-reads {output}
    """

rule generate_asvs_and_derep_tables:
    input:
        seqs = '{outdir}/simulation/asv_ref/{region}-seqs.qza',
        table = "{outdir}/simulation/samples/table.qza",
    output:
        joined_table="{outdir}/simulation/by_region/{region}/table.qza",
        joined_repset="{outdir}/simulation/by_region/{region}/rep_seqs.qza",
    params:
        length = lambda w: config['simulation'][w.region]['length'],
        seed = lambda w: config['simulation'][w.region]['asv_seed']
    shell:
        """
        python bin/simulate_tables_exact.py \
         --sequences {input.seqs} \
         --sample-table {input.table} \
         --read-length {params.length} \
         --random {params.seed} \
         --table {output.joined_table} \
         --reads {output.joined_repset}
        """


###############################################################################
#                           Set up the real dataset                           #
###############################################################################
#                                                                             #
# We're working with a subset of a set of benchmarking data from Hugerth et   #
# al. (DOI:L ). The original experiment compared the effect of multiple       #
# primer pairs and PCR techniques on taxonomic resolution in vaginal          #
# microbiome samples.                                                         #
#                                                                             #
# Here, we'll use the V1-3 primers from the human microbiome project (i hope) #
# without a chlamydia spikein, and V3-4 primers from Hugerth 2014 amplified   #
# with a 2 step PCR protocol.                                                 # 
#
# ...  
#                                                                             #
###############################################################################

###################################### v13 ####################################


rule trim_data_v13:
    input:
        config['real']['v13']['input_path']
    output:
        '{outdir}/real/by_region/v13/no-primers.qza'
    params:
        f_primer=config['real']['v13']['fwd']
    shell:
        """
        qiime cutadapt trim-single \
         --i-demultiplexed-sequences {input} \
         --p-front {params.f_primer} \
         --p-no-discard-untrimmed \
         --o-trimmed-sequences {output}
        """

rule denoise_v13:
    input:
        '{outdir}/real/by_region/v13/no-primers.qza'
    output:
        table='{outdir}/real/by_region/v13/pre-table.qza',
        seqs='{outdir}/real/by_region/v13/pre-rep_seqs.qza',
        stats='{outdir}/real/by_region/v13/stats.qza'
    shell:
        """
        qiime dada2 denoise-single \
         --i-demultiplexed-seqs {input} \
         --p-trunc-len 280 \
         --p-n-threads 4 \
         --o-table {output.table} \
         --o-representative-sequences {output.seqs} \
         --o-denoising-stats {output.stats}
        """

###################################### v34 ####################################

rule trim_data_v34:
    input:
        config['real']['v34']['input_path']
    output:
        '{outdir}/real/by_region/v34/no-primers.qza'
    params:
        f_primer=config['real']['v34']['fwd'],
        r_primer=config['real']['v34']['rev']
    shell:
        """
        qiime cutadapt trim-paired \
         --i-demultiplexed-sequences {input} \
         --p-front-f {params.f_primer} \
         --p-front-r {params.r_primer} \
         --p-no-discard-untrimmed \
         --o-trimmed-sequences {output}
        """

rule denoise_v34:
    input:
        '{outdir}/real/by_region/v34/no-primers.qza'
    output:
        table='{outdir}/real/by_region/v34/pre-table.qza',
        seqs='{outdir}/real/by_region/v34/pre-rep_seqs.qza',
        stats='{outdir}/real/by_region/v34/stats.qza'
    shell:
        """
        qiime dada2 denoise-paired \
         --i-demultiplexed-seqs {input} \
         --p-trunc-len-f 280 \
         --p-trunc-len-r 280 \
         --p-n-threads 4 \
         --o-table {output.table} \
         --o-representative-sequences {output.seqs} \
         --o-denoising-stats {output.stats}
        """



########################### Consistent trim length ############################

rule trim_posthoc:
    input:
        table='{outdir}/real/by_region/{region}/pre-table.qza',
        seqs='{outdir}/real/by_region/{region}/pre-rep_seqs.qza',
    output:
        table='{outdir}/real/by_region/{region}/table.qza',
        seqs='{outdir}/real/by_region/{region}/rep_seqs.qza',
    params:
        length = lambda w: config['real'][w.region]['length']
    shell:
        """
        qiime sidle trim-dada2-posthoc \
         --i-table {input.table}  \
         --i-representative-sequences {input.seqs} \
         --p-trim-length {params.length} \
         --p-hashed-feature-ids \
         --o-trimmed-table {output.table} \
         --o-trimmed-representative-sequences {output.seqs} \
        """

###############################################################################
#                         Format Database For Sidle                           #
###############################################################################
#                                                                             #
# ... 
#                                                                             #
###############################################################################

rule filter_silva_reference:
    input:
        seqs=config['simulation']['taxonomy']['rep_seq'],
        taxa=config['simulation']['taxonomy']['taxonomy']
    output:
        seqs=config['simulation']['taxonomy']['filt_rep_seq'],
    params:
        degen=config['simulation']['num_degen'],
        include=config['simulation']['include'],
        tmp='tmp_filt.qza'
    shell:
        """
        qiime sidle filter-degenerate-sequences \
         --i-sequences {input.seqs} \
         --p-max-degen {params.degen} \
         --p-n-workers 4 \
         --o-filtered-sequences {params.tmp}

        qiime taxa filter-seqs \
         --i-sequences {params.tmp} \
         --i-taxonomy {input.taxa} \
         --p-include {params.include} \
         --p-mode "contains" \
         --o-filtered-sequences {output.seqs} 

        rm {params.tmp}
        """


rule extract_sidle_regions:
    input:
        lambda w: config[w.dataset]['taxonomy']['filt_rep_seq']
    output:
        '{outdir}/{dataset}/database/{region}_extracted_reads.qza'
    params:
        fwd = lambda w: config[w.dataset][w.region]['fwd'],
        rev =lambda w: config[w.dataset][w.region]['rev'],
    shell:
        """
        qiime feature-classifier extract-reads \
         --i-sequences {input} \
         --p-f-primer {params.fwd} \
         --p-r-primer {params.rev} \
         --o-reads {output}
        """

rule prepare_extracted_region:
    input: 
        '{outdir}/{dataset}/database/{region}_extracted_reads.qza'
    output:
        kmers = '{outdir}/{dataset}/database/{region}_kmers.qza',
        map_ = '{outdir}/{dataset}/database/{region}_map.qza',
    params:
        fwd = lambda w: config[w.dataset][w.region]['fwd'],
        rev =lambda w: config[w.dataset][w.region]['rev'],
        trim = lambda w: config[w.dataset][w.region]['length'],
        region = lambda w: w.region
    shell:
        """
        qiime sidle prepare-extracted-region \
         --i-sequences {input} \
         --p-fwd-primer {params.fwd} \
         --p-rev-primer {params.rev} \
         --p-region {params.region} \
         --p-trim-length {params.trim} \
         --p-n-workers 2 \
         --o-collapsed-kmers {output.kmers} \
         --o-kmer-map {output.map_}
        """

###############################################################################
#                               Reference data                                #
###############################################################################
#                                                                             #
# ... 
#                                                                             #
###############################################################################

################################# Simulation ##################################

rule copy_sim_reference:
    input:
        table = '{outdir}/simulation/samples/table.qza', 
        taxa = config['simulation']['taxonomy']['reference']['taxonomy'],
        tree = config['simulation']['taxonomy']['reference']['tree'],
    output:
        table = '{outdir}/simulation/merged/reference/table.qza', 
        taxa = '{outdir}/simulation/merged/reference/taxonomy.qza',
        tree = '{outdir}/simulation/merged/reference/tree.qza',
    shell:
        """
        cp {input.table} {output.table}
        cp {input.taxa} {output.taxa}
        cp {input.tree} {output.tree}
        """

################################## Real Data ##################################

rule classify_regional_taxonomy:
    input:
        seqs='{outdir}/real/by_region/v34/pre-rep_seqs.qza',
        classifier=config['real']['taxonomy']['classifier']
    output:
        output='{outdir}/real/merged/reference/taxonomy.qza',
    shell:
        """
        qiime feature-classifier classify-sklearn \
         --i-reads {input.seqs} \
         --i-classifier {input.classifier} \
         --o-classification {output}
        """

rule get_real_reference_data:
    input:
        table='{outdir}/real/by_region/v34/pre-table.qza',
        taxonomy='{outdir}/real/merged/reference/taxonomy.qza',
    output:
        '{outdir}/real/merged/reference/table.qza',
    shell:
        """
        qiime taxa filter-table \
         --i-table {input.table} \
         --i-taxonomy {input.taxonomy} \
         --p-include ';c__' \
         --p-mode contains \
         --o-filtered-table {output}
        """

rule get_real_representative_seqs:
    input:
        table='{outdir}/real/merged/reference/table.qza',
        seqs='{outdir}/real/by_region/{region}/pre-rep_seqs.qza',
    output:
        '{outdir}/real/merged/reference/rep_seq.qza',
    shell:
        """
        qiime feature-table filter-seqs \
         --i-table {input.table} \
         --i-data {input.seqs} \
         --o-filtered-data {output}
        """

###############################################################################
#                                  Merge Data                                 #
###############################################################################
#                                                                             #
# ... 
#                                                                             #
###############################################################################

rule concatenate_single_asv_counts:
    input: 
        lambda w: config[w.dataset]['asv_merge']['tables']
    output:
        "{outdir}/{dataset}/merged/table.qza"
    run:
        from qiime2 import Artifact
        import qiime2.plugins.feature_table.actions as q2_table

        tables = [Artifact.load(fp) for fp in input]
        
        merged = q2_table.merge(tables, overlap_method='sum').merged_table
        merged.save(str(output))


rule concatenate_single_rep_seqs:
    input: 
        lambda w: config[w.dataset]['asv_merge']['sequences']
    output:
        "{outdir}/{dataset}/merged/rep_seqs.qza"
    run:
        from qiime2 import Artifact
        import qiime2.plugins.feature_table.actions as q2_table

        tables = [Artifact.load(fp) for fp in input]
        
        merged = q2_table.merge_seqs(tables).merged_data
        merged.save(str(output))


###############################################################################
#                                     OTUs                                    #
###############################################################################
#                                                                             #
# ... 
#                                                                             #
###############################################################################

rule cluster_otus:
    input:
        table="{outdir}/{dataset}/merged/table.qza",
        seqs="{outdir}/{dataset}/merged/rep_seqs.qza",
        ref=lambda w: config[w.dataset]['taxonomy']['filt_rep_seq'],
    output:
        table="{outdir}/{dataset}/merged/otus/table.qza",
        centroids="{outdir}/{dataset}/merged/otus/centroids-trimmed.qza",
        discard="{outdir}/{dataset}/merged/otus/discard-trimmed.qza",
    params:
        perc_identity=lambda w: config[w.dataset]['perc_identity'],
        threads=config['threads']
    shell:
        """
        qiime vsearch cluster-features-closed-reference \
         --i-sequences {input.seqs} \
         --i-table {input.table} \
         --i-reference-sequences {input.ref} \
         --p-perc-identity {params.perc_identity} \
         --p-threads {params.threads} \
         --o-clustered-table {output.table} \
         --o-clustered-sequences {output.centroids} \
         --o-unmatched-sequences {output.discard}
        """

rule get_otu_tree:
    input:
        lambda w: config[w.dataset]['taxonomy']['tree']
    output:
        "{outdir}/{dataset}/merged/otus/tree.qza"
    shell:
        """
        cp {input} {output}
        """

rule get_otu_taxonomy:
    input:
        lambda w: config[w.dataset]['taxonomy']['taxonomy']
    output:
        "{outdir}/{dataset}/merged/otus/taxonomy.qza"
    shell:
        """
        cp {input} {output}
        """


###############################################################################
#                                     ASVs                                    #
###############################################################################
#                                                                             #
# ... 
#                                                                             #
###############################################################################

rule classify_taxonomy:
    input:
        seqs="{outdir}/{dataset}/merged/rep_seqs.qza",
        classifier=lambda w: config[w.dataset]['taxonomy']['classifier']
    output:
        output='{outdir}/{dataset}/merged/asvs/taxonomy.qza',
    shell:
        """
        qiime feature-classifier classify-sklearn \
         --i-reads {input.seqs} \
         --i-classifier {input.classifier} \
         --o-classification {output}
        """

rule get_asv_table:
    input:
        table="{outdir}/{dataset}/merged/table.qza",
        taxonomy='{outdir}/{dataset}/merged/asvs/taxonomy.qza',
    output:
        '{outdir}/{dataset}/merged/asvs/table.qza',
    params:
        lambda w: config[w.dataset]['taxonomy']['class_str']
    shell:
        """
        qiime taxa filter-table \
         --i-table {input.table} \
         --i-taxonomy {input.taxonomy} \
         --p-include ";D_1__" \
         --p-mode contains \
         --o-filtered-table {output}
        """

rule get_asv_rep_seq:
    input:
        table='{outdir}/{dataset}/merged/asvs/table.qza',
        seqs="{outdir}/{dataset}/merged/rep_seqs.qza",
    output:
        '{outdir}/{dataset}/merged/asvs/rep_seq.qza',
    shell:
        """
        qiime feature-table filter-seqs \
         --i-table {input.table} \
         --i-data {input.seqs} \
         --o-filtered-data {output}
        """

rule build_insertion_tree:
    input:
        seqs="{outdir}/{dataset}/merged/rep_seqs.qza",
        reference=lambda w: config[w.dataset]['taxonomy']['sepp_reference']
    output:
        tree="{outdir}/{dataset}/merged/asvs/tree.qza",
        placements="{outdir}/{dataset}/merged/asvs/placements.qza"
    params:
        threads=config['threads']
    shell:
        """
        qiime fragment-insertion sepp \
            --i-representative-sequences {input.seqs} \
            --i-reference-database {input.reference} \
            --o-tree {output.tree} \
            --p-threads {params.threads} \
            --o-placements {output.placements}
        """

###############################################################################
#                                 Sidle                             #
###############################################################################
#                                                                             #
# The final approach is a kmer-based alignment (SMURF). The idea here is to   #
# calculate an alignment probability baesed on the relationhsip between a     #
# a sequence and the reference and then to solve the mixture over that region #
# to reconstruct a single table.                                              #
#                                                                             #
###############################################################################

rule align_sidle_regions:
    input:
        kmers='{outdir}/{dataset}/database/{region}_kmers.qza',
        rep_seq='{outdir}/{dataset}/by_region/{region}/rep_seqs.qza',
    output:
        alignment='{outdir}/{dataset}/by_region/{region}/sidle_alignment.qza',
    params:
        region = lambda w: w.region,
        mismatch = lambda w: config[w.dataset][w.region]['mismatch'],
        threads=config['threads']
    shell:
        """
        qiime sidle align-regional-kmers \
         --i-kmers {input.kmers} \
         --i-rep-seq {input.rep_seq} \
         --p-region {params.region} \
         --p-max-mismatch {params.mismatch} \
         --p-n-workers 4 \
         --p-chunk-size 75 \
         --o-regional-alignment {output.alignment} \
         --verbose
        """

rule reconstruct_table_single:
    input:
        kmers=lambda w: config[w.dataset]['asv_merge']['kmers'],
        align=lambda w: config[w.dataset]['asv_merge']['alignment'],
        tables=lambda w: config[w.dataset]['asv_merge']['tables'],
    output:
        table="{outdir}/{dataset}/merged/sidle/table.qza",
        summary="{outdir}/{dataset}/merged/sidle/reconstruction-summary.qza",
        db_map="{outdir}/{dataset}/merged/sidle/reconstruction-map.qza",
    params:
        lambda w: config[w.dataset]['regions'],
        threads=config['threads']
    run:
        from qiime2 import Artifact
        import qiime2.plugins.sidle.actions as q2_sidle

        num_regions = len(params[0])
        regions = params[0]
        nworkers = params[1]

        kmers = [Artifact.load(str(fp_)) 
                 for fp_ in input[:num_regions]]
        align = [Artifact.load(str(fp_)) 
                 for fp_ in input[num_regions:2*num_regions]]
        tables = [Artifact.load(str(fp_)) 
                    for fp_ in input[2*num_regions:]]

        table, summary, db_map = q2_sidle.reconstruct_counts(
            kmer_map=kmers,
            regional_alignment=align,
            regional_table=tables,
            region=regions,
            region_normalize='average',
            n_workers=3,
            # n_workers=nworkers,
            # verbose=True,
            )

        table.save(str(output[0]))
        summary.save(str(output[1]))
        db_map.save(str(output[2]))


rule reconstruct_taxonomy:
    input:
        db_map="{outdir}/{dataset}/merged/sidle/reconstruction-map.qza",
        taxa = lambda w: config[w.dataset]['taxonomy']['taxonomy'],
    output:
        "{outdir}/{dataset}/merged/sidle/taxonomy.qza",
    params:
        database = lambda w: config[w.dataset]['taxonomy']['database']
    shell:
        """
        qiime sidle reconstruct-taxonomy \
         --i-reconstruction-map {input.db_map} \
         --i-taxonomy {input.taxa} \
         --p-database {params.database} \
         --p-define-missing merge \
         --o-reconstructed-taxonomy {output}
        """

rule reconstruct_fragments_single:
    input:
        summary="{outdir}/{dataset}/merged/sidle/reconstruction-summary.qza",
        db_map="{outdir}/{dataset}/merged/sidle/reconstruction-map.qza",
        alignment = lambda w: config[w.dataset]['taxonomy']['aligned_repset'],
        kmers=lambda w: config[w.dataset]['asv_merge']['kmers'],
    output:
         "{outdir}/{dataset}/merged/sidle/reconstruction-fragments.qza",
    params:
        lambda w: config[w.dataset]['regions'],
    run:
        from qiime2 import Artifact
        import qiime2.plugins.sidle.actions as q2_sidle

        summary = Artifact.load(str(input[0]))
        db_map = Artifact.load(str(input[1]))
        aligned_seqs = Artifact.load(str(input[2]))
        kmer_maps = [Artifact.load(str(fp)) for fp in input[3:]]

        regions = params[0]

        representative_fragments = q2_sidle.reconstruct_fragment_rep_seqs(
            reconstruction_map=db_map,
            reconstruction_summary=summary,
            aligned_sequences=aligned_seqs,
            kmer_map=kmer_maps,
            region=regions
            ).representative_fragments
        representative_fragments.save(str(output))


rule reconstructed_fragment_insertion:
    input:
        fragments="{outdir}/{dataset}/merged/sidle/reconstruction-fragments.qza",
        reference = lambda w: config[w.dataset]['taxonomy']['sepp_reference']
    output:
        tree = "{outdir}/{dataset}/merged/sidle/tree.qza",
        placements = "{outdir}/{dataset}/merged/sidle/placements.qza"
    shell:
        """
        qiime fragment-insertion sepp \
            --i-representative-sequences {input.fragments} \
            --i-reference-database {input.reference} \
            --p-threads 6 \
            --o-tree {output.tree} \
            --o-placements {output.placements}
        """

###############################################################################
#                                    Diversity                                #
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################

rule rarify_table:
    input:
        table="{outdir}/{dataset}/merged/{method}/{table}.qza",
    output:
        table="{outdir}/{dataset}/merged/{method}/{table}-{depth}.qza",
    params:
        depth=lambda w: w.depth
    shell:
        """
        qiime feature-table rarefy \
         --i-table {input.table} \
         --p-sampling-depth {params.depth} \
         --o-rarefied-table {output.table}
        """

rule collapse_to_level:
    input:
        table="{outdir}/{dataset}/merged/{method}/table.qza",
        taxonomy="{outdir}/{dataset}/merged/{method}/taxonomy.qza",
    output:
        "{outdir}/{dataset}/merged/{method}/{level}-table.qza"
    params:
        collapse = lambda w: config['collapse'][w.level]
    shell:
        """
        qiime taxa collapse \
         --i-table {input.table} \
         --i-taxonomy {input.taxonomy} \
         --p-level {params} \
         --o-collapsed-table {output}
        """

rule rarefy_iteration:
    input:
        table="{outdir}/{dataset}/merged/{method}/table.qza",
    output:
        table="{outdir}/{dataset}/merged/{method}/rarified/{i}.qza",
    params:
        depth=10000
    shell:
        """
        qiime feature-table rarefy \
         --i-table {input.table} \
         --p-sampling-depth {params.depth} \
         --o-rarefied-table {output.table}
        """

# ################################# Alpha Diversity ###############################
        
rule alpha_iter:
    input:
        table="{outdir}/{dataset}/merged/{method}/rarified/{i}.qza",
        tree="{outdir}/{dataset}/merged/{method}/tree.qza",
    output:
        "{outdir}/{dataset}/merged/{method}/rarified-alpha/{metric}/{i}.qza",
    params:
        metric = lambda w: w.metric
    run:
        from qiime2 import Artifact
        import qiime2.plugins.diversity.actions as q2_diversity

        table = Artifact.load(str(input[0]))
        tree = Artifact.load(str(input[1]))
        metric = params[0]

        if metric == 'faith_pd':
            alpha = q2_diversity.alpha_phylogenetic(
                table=table,
                phylogeny=tree,
                metric=metric,
                ).alpha_diversity
        else:
            alpha = q2_diversity.alpha(
                table=table,
                metric=metric,
                ).alpha_diversity
        alpha.save(output[0])

# ################################# Beta Diversity ###############################

rule beta_phylogenetic:
    input:
        table="{outdir}/{dataset}/merged/{method}/rarified/{i}.qza",
        tree="{outdir}/{dataset}/merged/{method}/tree.qza",
    output:
        "{outdir}/{dataset}/merged/{method}/rarified-beta/{weight}-unifrac/{i}.qza",
    params:
        metric=lambda w: '{}_unifrac'.format(w.weight)
    shell:
        """
        qiime diversity beta-phylogenetic \
         --i-table {input.table} \
         --i-phylogeny {input.tree} \
         --p-metric {params.metric} \
         --o-distance-matrix {output}
        """

rule genus_braycurtis:
    input:
        table="{outdir}/{dataset}/merged/{method}/rarified/{i}.qza",
        taxonomy="{outdir}/{dataset}/merged/{method}/taxonomy.qza",
    output:
        dm="{outdir}/{dataset}/merged/{method}/rarified-beta/genus-braycurtis/{i}.qza",
        table="{outdir}/{dataset}/merged/{method}/rarified-genus/{i}.qza",
    shell:
        """
        qiime taxa collapse \
         --i-table {input.table} \
         --i-taxonomy {input.taxonomy} \
         --p-level 6 \
         --o-collapsed-table {output.table}

        qiime diversity beta \
        --i-table {output.table} \
        --p-metric braycurtis \
        --o-distance-matrix {output.dm}
        """

rule braycurtis:
    input:
        table="{outdir}/{dataset}/merged/{method}/rarified/{i}.qza",

    output:
        "{outdir}/{dataset}/merged/{method}/rarified-beta/braycurtis/{i}.qza",
    shell:
        """
        qiime diversity beta \
        --i-table {input.table} \
        --p-metric braycurtis \
        --o-distance-matrix {output}
        """


rule adonis:
    input:
        dm = "{outdir}/{dataset}/merged/{method}/rarified-beta/{metric}/{i}.qza",
        meta='{outdir}/simulation/samples/metadata.tsv',
    output:
        '{outdir}/{dataset}/merged/{method}/rarified-beta/{metric}-adonis/{i}.tsv'
    params:
        dir_="{outdir}/{dataset}/merged/{method}/rarified-beta/{metric}/{i}/dm"
    shell:
        """
        qiime tools export \
         --input-path {input.dm} \
         --output-path {params.dir_}

        Rscript bin/adonis_single.R \
         -d {params.dir_}/distance-matrix.tsv \
         -m {input.meta} \
         -o {output}

        rm -r {params.dir_}
        """


###############################################################################
#                                Results                                      #
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################

rule alpha_table:
    input:
        expand("{outdir}/simulation/merged/{method}/rarified-alpha/{metric}/{i}.qza",
               outdir=config['output_dir'],
               method=['reference', 'otus', 'asvs', 'sidle'],
               metric=['observed_features', 'shannon', 'pielou_e', 'faith_pd'],
               i=[x for x in np.arange(0, config['rare_iter'])]),
        'data/output/simulation/samples/metadata.tsv',
    output:
        "data/output/tables/table-1-alpha.tsv"
    notebook:
        "ipynb/Table1-CompareAlpha.ipynb"

rule beta_table:
    input:
        expand("{outdir}/{dataset}/merged/{method}/rarified-beta/{metric}/{i}.qza",
               outdir='ipynb/data/output',
               dataset=['simulation'],
               method=['reference', 'otus', 'asvs', 'sidle'],
               metric=['braycurtis', 'genus-braycurtis', 'unweighted-unifrac', 'weighted-unifrac'],
               i=[x for x in np.arange(0, config['rare_iter'])]
               ),
        expand('{outdir}/{dataset}/merged/{method}/rarified-beta/{metric}-adonis/{i}.tsv',
               outdir='ipynb/data/output',
               dataset=['simulation'],
               method=['reference', 'otus', 'asvs', 'sidle'],
               metric=['braycurtis', 'genus-braycurtis', 'unweighted-unifrac', 'weighted-unifrac'],
               i=[x for x in np.arange(0, config['rare_iter'])]
               ),
        'data/output/simulation/samples/metadata.tsv',
    output:
        "data/output/tables/table_2_beta_rarefaction.tsv"
    notebook:
        "Table2-CompareBeta.ipynb"

rule compare_taxonomy:
    input:
        'data/output/simulation/samples/metadata.tsv',
        expand('{outdir}/{dataset}/merged/{method}/{artifact}.qza',
               **config['output']['tables']
               ),
    output:
        'data/output/figures/figure2/fig2_class_relative_abund.png',
        'data/output/tables/table_s2_class_comparison.tsv'
    notebook:
        'ipynb/Figure1-Taxonomy.ipynb'

rule real_data_figure:
    input:
        expand('{outdir}/{dataset}/merged/{method}/{artifact}.qza',
               **config['output']['tables']
               ),
        'data/inputs/real/manifest-v13.tsv',
    output:
        'data/output/figures/figure3/figure2_vaginal.png'
    notebook:
        'ipynb/Figure2-RealData.ipynb'
