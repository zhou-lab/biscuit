.. BISCUIT documentation master file, created by
   sphinx-quickstart on Sun Apr 10 10:26:36 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BISCUIT User Guide
===========================

BISCUIT is a software suite for analyzing sodium bisulfite conversion-based DNA methylation data.

Download and Install
######################

The latest release is available `here <https://github.com/zwdzwd/biscuit/releases/latest>`_.

To install, simply unzip and make,

.. code:: bash
	  
   unzip download
   make

The binaries will be created at bin/.

Contents:

.. toctree::
   :maxdepth: 2

   align_bisulfite_converted_reads
   measure_cytosine_retention_and_snp
   duplicate_marking
   nucleosome_depleted_region


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

