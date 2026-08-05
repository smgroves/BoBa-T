============================================================================
BoBa-T: A single-cell gene regulatory network inference and simulation tool
============================================================================

.. image:: assets/bobaT_logo.png
   :align: left
   :width: 200px
   :alt: BoBa-T logo

**BoBa-T** (**Bo**olean **Ba**yesian **T**ranscription Factor Networks) is a suite of network inference tools to derive and simulate gene
regulatory networks from transcriptomics (and chromatin accessibility) data. It is the expanded and refined single-cell update to
BooleaBayes, published in PLOS Computational Biology, `Wooten et al. 2019`_. The benefit of using BoBa-T as opposed to BooleaBayes is that it can infer regulatory networks from **single-cell RNA-seq data**, and it can incorporate **chromatin accessibility data (scATAC-seq)** to improve the accuracy of the inferred network. The method is designed to be robust to noise in single-cell data, and can handle large datasets with many genes and cells. See the paper here: Bhattacharya et al. (2026)XXX. 


Why use BoBa-T?
==================

A few benefits of using BoBa-T compared to other GRN methods:

* BoBa-T infers **a single network that can describe all stable states of a system**, rather than different networks for each cell state (see the underlying theory, based in dynamical systems theory, in papers such as `Huang et al. (2009) <https://pmc.ncbi.nlm.nih.gov/articles/PMC2754594/>`_). This allows you to perturb the single network *in silico* to make predictions about cell state transitions. 
* **Because the rules are Boolean^, they are interpretable and testable.** Each edge in the network represents a direct, physical regulatory interaction between a transcription factor and its target gene. Nonlinear relationships between transcription factors and their targets can be captured in the inferred rules, but only direct TF-target relationships are represented in the network.
* Similarly, because BoBa-T finds a Boolean network that describes the dynamics of the system, it can be used to identify *attractors (stable states) of the system*, and map observed cells onto these attractors. This allows for the identification of cell types and states, and the prediction of how cells may transition between these states. Therefore, **the steady states of the system are inferred from the underlying regulatory logic directly,** whereas many other methods infer steady states from the data itself, which can be confounded by noise and batch effects. Importantly, this means that **steady states can be predicted even when they are not observed in the data**, which is particularly useful for single-cell data, where some cell types may be rare or missing from the dataset.
* Many GRN inference methods find correlations between genes, but do not infer the underlying regulatory logic. BoBa-T infers the regulatory logic of each transcription factor as a pseudo-Boolean function of its parent nodes, which allows for the **inference of complex regulatory logic, such as AND, OR, and NOT relationships between transcription factors and their target genes.** This is particularly useful as an interpretable and testable model of gene regulation. Other GRN inference methods that find regulatory relationships often find a single edge weight between a TF and its target that applies to all conditions; **BoBa-T inherently infers interactions between TFs**. For example, if both TF1 and TF2 regulate GeneX, BoBa-T can infer that TF1 only regulates GeneX when TF2 is also present, or that TF1 only regulates GeneX when TF2 is absent. This allows for the inference of combinatorial regulation, which is a key feature of gene regulatory networks (`Buchler et al. (2003) <https://doi.org/10.1073/pnas.0930314100>`_ and `Balaji et al. (2006) <https://doi.org/10.1016/j.jmb.2006.04.029>`_ ).
* While BoBa-T is Boolean in nature, Boolean rules are often underfit due to the sparsity of single-cell data. Compared to other Boolean network inference methods, **BoBa-T is design to be robust to this sparsity using Bayesian principles to infer the most likely regulatory logic given the data.** This allows for the inference of regulatory logic even when the data is sparse or noisy, which is often the case in single-cell RNA-seq data.
* To see the benefit of these choices in BoBa-T, we have benchmarked the method against other GRN inference methods on both simulated and real data at multiple comparison points. See Bhattacharya et al. (2026)XXX for details.


Package structure
======================

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


.. _Wooten et al. 2019: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007343

^The network rule inference is based on Boolean network theory, where the regulatory logic of each transcription factor is inferred as a pseudo-Boolean function of its parent nodes. This allows for the inference of complex regulatory logic, such as AND, OR, and NOT relationships between transcription factors and their target genes. This is particularly useful as an interpretable and testable model of gene regulation; each edge in the network represents a direct, physical regulatory interaction between a transcription factor and its target gene. Nonlinear relationships between transcription factors and their targets can be captured in the inferred rules, but only direct TF-target relationships are represented in the network. For an introduction to the underlying Boolean theory of gene regulation, see our original BooleaBayes paper `Wooten et al. 2019`_ and the following sources: `Introduction to Systems Biology by Uri Alon <https://www.amazon.com/dp/1439837171?lv=shuf&channelId=500&plpRedirect=mhFallback>`_ and `Bornholdt (2008) <https://doi.org/10.1098/rsif.2008.0132.focus>`_.
