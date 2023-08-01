.. _MCOMP:

Micro-C Comparative Analyses
============================

Introduction
------------

Biological questions are seldom answered by analysing single samples in isolation. It is often the case that an experiment aims to make comparisons between two (or more) biological conditions,
such as:

1)	Untreated wild type vs treatment
2)  Wild type vs knockout
3)  Normal sample vs tumor

In all cases the goal is to produce a list of differentially interacting regions in one condition relative to the other. The main output for comparative analses is analogous to what is expeected for differential gene expression,
where the primary result is a table of regions, the fold change between conditions, and a statistical measure of signficance. For Micro-C, we aim to identify regions of differential interaction directly from the matrix files. See
:ref:`previous steps <GCM>` to generate the required matrices for differential analysis.

Figure 1:

.. image:: /images/CA_MC_fig1.png


Differential Analysis
---------------------

**Question:** How do I perform differential analyses for Micro-C experiments?

**Process:** Mcool files are first converted to text files of a perferred resolution, and then used as input to the HiCcompare algorithm.

**Results:** Final results consist of a table of differentially interacting regions, fold change, and measure of statistical signficance.

**Files and tools needed:**
  - .cool, .mcool, .hic, or Hic-Pro files for each replicate and sample condition
  - `HiCcompare <https://www.bioconductor.org/packages/release/bioc/html/HiCcompare.html>`_ for single-replicate analysis or `multiHiCcompare <https://www.bioconductor.org/packages/release/bioc/html/multiHiCcompare.html>`_ for multiple replicate experiments.

As the design of differential analysis experiments are unique to each biological question, there are multiple possibilites for how the analysis can be set up. A common scenario is to compare two conditions
where each condition has two replicates, and is described in the `multiHiCcompare vignette <https://bioconductor.org/packages/devel/bioc/vignettes/multiHiCcompare/inst/doc/multiHiCcompare.html>`_. The HiCcompare package also contains
functions for conversion of various input files


**Interpreting results:**

Micro-C differential analysis produces a number of intermediate files in addition to the final results table. There are two main outputs to consider:

1) MD normalization plots
2) Differential regions table

MD is a concept introduced by the HiCcompare developers and is analogous to the Tukey's mean/difference plot. M corresponds to the log2 fold change between the two conditions, and D is the distance between the two interacting regions.
Loess normalization aims to eliminate the bias introduced by the influence of interaction distance on fold change bewteen two conditions. It is often useful to visualize the effect of normalization between conditions to ensure
the data is appropriate for downstream difference detection. An example effect of normalization is given below:

Figure 2:

.. image:: /images/CA_MC_fig2.png

For difference deteciton, the resulting output file is highly similar to what is expected for gene expression studies, where regions are listed and prioritized by a combination of fold change and a measure of statistical signficance.
Below is an example output from HicCompare:

.. csv-table::
   :file: tables/CA_MC_results.csv
   :header-rows: 1
   :widths: 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12
   :class: tight-table

The most relevant fields from the output will be:
  - adj.M -- the log fold change in coverage between the two conditions
  - p.adj -- a p-value, after correction for multiple hypothesis testing, on the statistical signficance of the observed fold change

**Considerations:**

   - Replication – It is generally advisable to have technical replicates for differential analyses, as this will produce more statistically robust results.
