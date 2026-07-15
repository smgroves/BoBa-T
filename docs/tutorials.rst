=========
Tutorials
=========

These tutorials walk through the full BoBa-T workflow, from building a base
network structure to running *in silico* perturbations. Each tutorial contains
instructions to run the example data, the expected output, and explanations of
the code.

The workflow has four main stages, followed by a set of validation and
downstream-analysis tools:

#. **Build a base network structure** — define which transcription factors can
   regulate which, either from scATAC-seq data or from the original BooleaBayes
   ChEA-based approach.
#. **Infer network rules** — learn the Boolean regulatory logic for each node,
   from either two timepoints or pseudotime.
#. **Find attractors** — identify the stable states (phenotypes) of the network.
#. **Run in silico perturbations** — predict how knockouts and overexpression
   reshape the attractor landscape.

.. contents:: On this page
   :local:
   :depth: 1


1. Building a base network structure
=====================================

The base network defines the candidate regulatory edges (which transcription
factors may regulate which genes) before rule inference. BoBa-T supports two
ways to build it:

* **From scATAC-seq data** — derive edges from chromatin accessibility using a
  DIRECT-NET input, then prune with LASSO feature selection.
* **From ChEA databases (original BooleaBayes method)** — build the network from
  curated transcription-factor–target interactions.

.. toctree::
   :maxdepth: 1

   tutorials/network_example
   tutorials/network_chea


2. Inferring network rules
==========================

Rule inference (BoBa-T rule generation) learns the Boolean update logic for each
node in the network from binarized expression data. Rules can be fit from either
of two experimental designs:

* **Two timepoints** — pair each cell/state with a later timepoint.
* **Pseudotime** — use an RNA-velocity / pseudotime ordering to define the
  direction of change.

.. toctree::
   :maxdepth: 1

   tutorials/inference_example
   tutorials/inference_two_timepoints
   tutorials/inference_pseudotime


3. Finding attractors
=====================

Once rules are fit, BoBa-T identifies the network's attractors — the stable
states that correspond to cell phenotypes — and maps observed cells onto them.

.. toctree::
   :maxdepth: 1

   tutorials/attractors


4. In silico perturbations
==========================

With rules and attractors in hand, BoBa-T predicts how genetic perturbations
(knockouts and overexpression) destabilize or shift the attractor landscape.

.. toctree::
   :maxdepth: 1

   tutorials/perturbations


Validation and downstream analysis
===================================

Peripheral tools for evaluating rule accuracy and exploring the dynamics of the
fitted network, including random-walk trajectories.

.. toctree::
   :maxdepth: 1

   tutorials/validation
   tutorials/trajectories
