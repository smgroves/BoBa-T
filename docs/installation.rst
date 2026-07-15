============
Installation
============

System requirements
====================

BoBa-T is compatible with Python 3.8 and above and runs on CPU hardware. It has
been tested on Windows and macOS operating systems.

Dependencies
============

The ``graph-tool`` Python package must be installed separately — it cannot be
installed with ``pip``. Install it with ``conda`` or ``homebrew`` as described in
the `graph-tool installation instructions
<https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions>`_.

All other dependencies are installed automatically with the package and are
listed in ``setup.py``.

Installing BoBa-T
=================

Install the latest release from PyPI::

    pip install bobaT

Typical install time is less than one minute on a standard laptop.

To install from source::

    git clone https://github.com/smgroves/BoBa-T.git
    cd BoBa-T
    pip install -e .
