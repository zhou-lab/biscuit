****************************************
Detect Nucleosome Depleted Region
****************************************

The `ndr` subcommand in BISCUIT detects nucleosome depleted region (NDR) from Nucleosome Occupancy and Methylome Sequencing (NOMe-seq) data. To call nucleosome depleted region, one issues

.. code:: bash

   biscuit vcf2bed -t gch -c HCT116.vcf.gz | biscuit-develop ndr -b - -o OUT

which infer on each base of GCH context a state whether the region is depleted of occupancy. 

::

   ...
   chr18   11566   1
   chr18   11574   1
   chr18   11575   1
   chr18   11687   1
   chr18   11704   0
   chr18   11705   0
   ...

`-c` option further collapse the state calls into regions:

.. code:: bash

   biscuit vcf2bed -t gch -c HCT116.vcf.gz | biscuit-develop ndr -b - -o OUT

outputs regions spanned by consecutive depletion regions.

::

   ...
   chr18   11515   11516
   chr18   11524   11687
   chr18   11902   11945
   ...
