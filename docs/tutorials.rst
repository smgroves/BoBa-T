=========
Tutorials
=========

These tutorials walk through the full BoBa-T workflow, from building a base network structure to running *in silico* perturbations. Each tutorial contains instructions to run the example data, the expected output, and explanations of the code.

The workflow has four main stages, as well as a set of validation, plotting, and downstream-analysis tools. 

#. **Build a base network structure** — define which transcription factors can regulate others, either from scATAC-seq data for your given system of interest or using the original BooleaBayes approach that pulls information from ChIP-seq databases.
#. **Infer network rules** — learn the Boolean regulatory logic for each node. 
#. **Find attractors** — identify the stable states (phenotypes) of the network. 
#. **Run in silico perturbations** — predict how knockouts and overexpression
   reshape the attractor landscape.

.. contents:: On this page
   :local:
   :depth: 1


0. Input data
=====================================

* **scRNA-seq input** - BoBa-T expects scRNA-seq data. 

   * You should input a csv file of a samples (row) X genes (columns) matrix. 
   * For an example pipeline used in Bhattacharya et al., see this tutorial: Preprocessing scRNA-seq DataXXX. However, it is important to note that **BoBa-T is NOT a scRNA-seq preprocessing pipeline.** 
   * These normalized RNA counts are used optionally in step 1 (for LASSO regression) and in step 2: inferring network rules. 

* **scATAC-seq input** - To take full advantage of the method, BoBa-T can use scATAC-seq data for base network generation. 

   * The network structure generation (step 1) expects putative link information as a csv where each row is ``TF regulator, TF target, evidence[optional]``. 
   * An example of generation of this csv can be found in this tutorial: Preprocessing scATAC-seq Data XXX. However, it is important to note that **BoBa-T is NOT a scATAC-seq preprocessing pipeline.** 
   * BoBa-T was tested on multiomic data where the same cells have chromatin accessibility and RNA count data. While this is not required, we have not tested the algorithm on independent data sourcesXXX. 
   * In the absence of chromatin accessibility data, the method can fall back to the ChIP-seq database-searching method described in Wooten, Groves, et al (2019). A literature-derived putative gene list is useful for removing unrelated TFs, though not required.

.. list-table::
   :header-rows: 1

   * - scATAC-seq available?
     - scRNA-seq available?
     - pseudotime available?
     - multiple timepoints?
     - Network Generation Method
     - Rule-Fitting Method
   * - ❌ 
     - ❌
     - ❌
     - ❌
     - Use our previous `BooleaBayes <https://github.com/smgroves/BooleaBayes>`_ package to generate network from ChIP-seq databases
     - Use our previous `BooleaBayes <https://github.com/smgroves/BooleaBayes>`_ package to fit rules from bulk RNA-seq data
   * - ❌ 
     - ✅
     - ❌
     - ❌
     - Use the :doc:`base network generation method <tutorials/network_chea>`
     - Use the :doc:`base rule-fitting method <tutorials/inference_example>`
   * - ✅ 
     - ✅
     - ❌
     - ❌
     - Use the :doc:`chromatin accessibility and LASSO pruning method <tutorials/network_example>`
     - Use the :doc:`base rule-fitting method <tutorials/inference_example>`   
   * - ✅ 
     - ✅
     - ✅
     - ❌
     - Use the :doc:`chromatin accessibility and LASSO pruning method <tutorials/network_example>`
     - Use the :doc:`rule-fitting method <tutorials/inference_example>` 
   * - ✅ 
     - ✅
     - ❌
     - ✅
     - Use the :doc:`chromatin accessibility and LASSO pruning method <tutorials/network_example>`
     - Use the :doc:`pseudotime rule-fitting method <tutorials/inference_pseudotime>` 
   * - ✅ 
     - ✅
     - ✅
     - ❌
     - Use the :doc:`chromatin accessibility and LASSO pruning method <tutorials/network_example>`
     - Use the :doc:`two-timepoint rule-fitting method <tutorials/inference_two_timepoints>` 


1. Building a base network structure
=====================================

The base network defines the candidate regulatory edges (which transcription factors may regulate which genes) before rule inference. BoBa-T supports two ways to build it:

* **From scATAC-seq data** — derive edges from chromatin accessibility using a putative TF-target input, then prune with LASSO feature selection. This option uses scRNA-seq data for LASSO regression.

* **From ChEA databases (adaptation of the original BooleaBayes method)** — build the network from curated transcription-factor–target interactions.

There are options to preserve non-TF genes as readouts of the network for later interpretability; these genes will not affect the dynamics of the network.

Please note that unlike many other common GRN inference methods, BoBa-T infers pseudo-Boolean rules for each TF. This means that it is inferring ~2^N values for each TF, where N is the number of parent nodes in the network. Therefore, it is important to prune spurious edges *before* inferring network rules step 2. In practice, less then 10 parent nodes per TF will run in a reasonable amount of time. If step 2 below seems to be computationally inefficient, you should revisit this step for advancing pruning. Methods to do so using LASSO regression are shown in the first tutorial below.

.. toctree::
   :maxdepth: 1

   tutorials/network_example
   tutorials/network_chea


2. Inferring network rules
==========================

Rule inference (BoBa-T rule generation) learns the pseudo-Boolean update logic for each
node in the network from binarized expression data. Rules can be fit from either
of two experimental designs:

* **Single timepoint** — the rules are fit such that the current target expression is predicted from the current TF parent node expression profile. This makes the assumption that a large proportion of the dataset is already at steady state, and is useful for identifying rules of transitions when transition states cannot be captured.
* **Two timepoints** — pair each cell/state with a later timepoint. 

   * **Pseudotime** — use an RNA-velocity / pseudotime ordering to define the direction of change.

The default is to use a single timepoint of scRNA-seq data as input. You may also use data from two timepoints or pseudotime (e.g. run RNA velocity analyses first and then build the network on the true data + an inferred timepoint). Using two timepoints can give a sense of scaling to the temporal axis of the perturbation analysis, but a large difference in time scale between the timepoints and the time scale of TF regulation can lead to errors. 

.. toctree::
   :maxdepth: 1

   tutorials/inference_example
   tutorials/inference_two_timepoints
   tutorials/inference_pseudotime


3. Finding attractors
=====================

Once rules are fit, BoBa-T identifies the network's attractors — the stable
states that correspond to cell phenotypes — and maps observed cells onto them.

BoBa-T was designed to build *a single network* that can describe all stable states of a system, rather than different networks for each cell state. This allows you to perturb the single network to make predictions about cell state transitions, but also means that it implicitly assumes that cell states share a genetic background and cells can transition across the landscape. For example, you should not include data from multiple cell types that cannot reasonably transition from one state to another, such as including cancer cell populations with immune cell populations.

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
