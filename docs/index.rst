============================================================================
BoBa-T: A single-cell gene regulatory network inference and simulation tool
============================================================================

**BoBa-T** is a suite of network inference tools to derive and simulate gene
regulatory networks from transcriptomics data. It is the expanded and refined single-cell update to
BooleaBayes, published in PLOS Computational Biology,
`Wooten, Groves et al. (2019) <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007343>`_.

The benefit of using BoBa-T is that it can infer regulatory networks from **single-cell RNA-seq data**, and it can incorporate **chromatin accessibility data (scATAC-seq)** to improve the accuracy of the inferred network. The method is designed to be robust to noise in single-cell data, and can handle large datasets with many genes and cells. See the paper here: Bhattacharya et al. (2026)XXX. 

Importantly, BoBa-T infers *a single network* that can describe **all stable states of a system**, rather than different networks for each cell state (see the underlying theory, based in dynamical systems theory, in papers such as `Huang et al. (2009) <https://pmc.ncbi.nlm.nih.gov/articles/PMC2754594/>`_). This allows you to perturb the single network *in silico* to make predictions about cell state transitions. 

The network rule inference is based on Boolean network theory, where the regulatory logic of each transcription factor is inferred as a pseudo-Boolean function of its parent nodes. This allows for the inference of complex regulatory logic, such as AND, OR, and NOT relationships between transcription factors and their target genes. This is particularly useful as an interpretable and testable model of gene regulation; each edge in the network represents a direct, physical regulatory interaction between a transcription factor and its target gene. Nonlinear relationships between transcription factors and their targets can be captured in the inferred rules, but only direct TF-target relationships are represented in the network. For an introduction to the underlying Boolean theory of gene regulation, see our original BooleaBayes paper (Wooten, Groves et al. 2019) and the following sources: `Introduction to Systems Biology by Uri Alon <https://www.amazon.com/dp/1439837171?lv=shuf&channelId=500&plpRedirect=mhFallback>`_ and `Bornholdt (2008) <https://doi.org/10.1098/rsif.2008.0132.focus>`_.

The package is organized into the following modules:

* ``net`` — make or modify network structure
* ``load`` — loading data
* ``proc`` — processing
* ``rw`` — random walk
* ``plot`` — plotting
* ``tl`` — tools
* ``utils`` — utilities

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
