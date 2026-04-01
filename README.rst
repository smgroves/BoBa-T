=======================================================
BooleaBayes
=======================================================
.. image:: https://badge.fury.io/py/booleabayes.svg
    :target: https://pypi.org/project/booleabayes/
    :alt: Latest PYPi Version

BoBa-T is a suite of network inference tools to derive and simulate gene regulatory networks from transcriptomics data; it is our single-cell update to BooleaBayes, which was published in PLOS Computational Biology, `Wooten, Groves et al. (2019) <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007343>`_. 

System Requirements
~~~~~~~~~~~~~~~~~~~~~~~~

BoBa-T is compatible with Python 3.8 and above and runs on CPU hardware. It has been tested on Windows and macOS operating systems.

Dependencies
---------------------

The ``graph-tool`` python package will need to be installed. This can be installed with `Conda`, `homebrew`, etc as seen `here <https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions>`_. 

All other dependencies will be installed with this package and can be found within the `setup.py` file. 


Installation Guide
~~~~~~~~~~~~~~~~~~~~~~~~

To install ``boba-T``, please use::

    pip install bobaT

Typical install time is less than 1 minute on a standard laptop.


Instructions for Use:
~~~~~~~~~~~~~~~~~~~~~~~~

The BoBa-T package is organized into the following modules:

* ``net`` = make or modify network structure
* ``load`` = loading data
* ``proc`` = processing
* ``rw`` = random walk
* ``plot`` = plotting
* ``tl`` = tools
* ``utils`` = utilities

Demo:
~~~~~~~~~~~~~~~~~~~~~~~~

For more details on how to use these functions, please see the tutorials for network construction (`network_example.ipynb`) and inference (`inference_example.ipynb`). These tutorials contain instructions to run the data, expected output, and explanations of the code.
The network generation tutorial should run in less than 10 minutes on a standard laptop.
The inference tutorial should run in less than 30 minutes on a standard laptop.


See this `repository <https://github.com/smgroves/Bhattacharya2026>`_  for examples and documentation of how to use this package.