============
Key Concepts
============

This page defines, in plain language, the core vocabulary used throughout the tutorials, plots, and API — in particular the terms **attractor**, **basin**, **perturbation**, and **destabilization score**. It's meant as a quick working reference, not a substitute for the full theoretical treatment in `Wooten et al. 2019 <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007343>`_ (the original BooleaBayes paper) and Bhattacharya et al. (under review, the BoBa-T paper).

.. contents:: On this page
   :local:
   :depth: 1

Node, rule, and state
======================

* A **node** is a gene (usually a transcription factor) in the network. 
* A **rule** is the piece of regulatory logic BoBa-T infers for one node: given the current on/off (or continuous, 0-1) values of that node's regulators, the rule predicts that node's own value. Rules can encode AND/OR/NOT-style combinatorial logic between regulators, not just independent effects (see the Boolean-rule footnote in :doc:`index`).
* A **state** is a snapshot of every node's value at once — the full on/off profile of the whole network for one cell (or one point along a simulated walk).

Attractors and basins
======================

* An **attractor** is a *stable* state: once the network's dynamics reach it, they stay there. Concretely, every node's rule agrees with the node's own current value, so nothing has a reason to flip. Attractors are what BoBa-T maps onto observed cell phenotypes (e.g. "NE", "nonNE", "Intermediate").
* Because the network can be simulated from *any* starting state, not just the ones actually observed in your data, BoBa-T can find attractors that correspond to rare or entirely missing cell states — this is the sense in which "steady states can be predicted even when they are not observed in the data" (see :doc:`index`).
* The **basin of attraction** for a given attractor is the surrounding region of state-space that flows into it: if you start a random walk anywhere inside the basin and let it run with no perturbation, it will (with high probability) settle into that attractor. Basin membership is what :doc:`tutorials/attractors` uses to filter and label attractors, and what a "radius" parameter (e.g. ``tf_basin``) controls.

  .. note::
     Because BooleaBayes and BoBa-T are Bayesian-adapted Boolean network, the term **"pseudo-attractor"** is used to describe a state where the probability of leaving it *in any direction* is lower than the probability of staying. In other words, there is a greater than 0.5 probability of remaining in the attractor rather than flipping a node and moving away.

Perturbations
======================

A **perturbation** is an *in silico* intervention on the network: force one node's value to be fixed OFF (a **knockdown**) or fixed ON (an **activation**/overexpression) for the whole simulation, instead of letting its rule update it normally. Everything else in the network still follows its own rule as usual. This lets you ask "if I knocked out or overexpressed this one gene, where would the cell end up?" — including whether it leaves its current attractor's basin altogether and, if so, which attractor it moves toward instead.

Destabilization score
======================

The **destabilization score** for a gene, with respect to a given attractor, measures how much perturbing that gene shortens the time a random walk takes to *leave* that attractor's basin, relative to the unperturbed (control) case:

.. code-block:: text

   destabilization score = (mean steps-to-leave, perturbed − mean steps-to-leave, control) / mean steps-to-leave, control

A **negative** score means the perturbation makes cells leave the basin *faster* than normal — i.e. the perturbed gene is *destabilizing* that attractor/state. A score near zero means the perturbation has little effect on how long the cell stays in that state.

.. important::
   **Naming to watch for**: BoBa-T uses ``"Regulators"`` for genes whose **knockdown** significantly destabilizes a state, and ``"Destabilizers"`` for genes whose **activation** significantly destabilizes it. 

A note on "activator" / "repressor"
======================================

It's tempting to relabel destabilizing genes as "activators" or "repressors" of a state, but those terms already have an established, *different* meaning in TF biology: whether a transcription factor turns a specific target gene's expression up or down. A destabilization score describes something else — its effect on a *whole-network attractor*, not on any single target gene — and a gene can be a strong destabilizer of a state regardless of whether it acts as a transcriptional activator or repressor on its individual targets. Conflating the two would misrepresent what the score measures, so this documentation uses "destabilizing" / "stabilizing" (with respect to a named attractor) rather than "activator" / "repressor" throughout.
