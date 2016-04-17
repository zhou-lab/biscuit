**************************************************
Mark PCR Duplicates in Bisulfite-sequencing
**************************************************

Compared with non-bisulfite treated sequencing, PCR-duplication can be further resolved using bisulfite conversion strand information (T-rich or A-rich) since sodium-bisulfite conversion occurs before PCR amplification. In BISCUIT, one read is not considered the duplicate of another if the two reads are on different conversion strands. 

The following command mark duplicates in before.bam

.. code:: bash

   biscuit markdup before.bam after.bam

BISCUIT considers two paired-end reads as duplication if 

  + the mates are on the same chromosome;
  + the insert size is smaller than a threshold (defined by `-l` option);
  + Both mates are mapped to the same location;
  + Both inserts are with the same conversion strand.

All reads/inserts are marked as duplication except the one with the highest base quality sum. `-u` option let mate-unmapped paired-end reads be treated as pure single-end reads. By default, BISCUIT does not consider mate-unmapped paired-end reads as purely single-end reads. This further utilizes the read position in pair information in resolving PCR duplication (when two reads are both mate-unmapped, first read in pair is not the duplicate of the second read in another pair even if the coordinates match). Dangling reads are treated as single-end reads. Secondary mapping and unmapped reads are not processed and are directly flushed to output. `-r` option allows duplicate reads be removed instead of duplicate-marked. 

Supported tags
^^^^^^^^^^^^^^^^^^^^^

Conversion strand information must be encoded as a tag in the bam file. BISCUIT supports bisulfite conversion strand information inferred from multiple existing aligners, including

+----------+-----------+-----------------+
| aligner  | tag name  | tag candidates  |
+----------+-----------+-----------------+
| BWA-meth | YD        | f/r             |
+----------+-----------+-----------------+
| Bismark  | XG        | CT/GA           |
+----------+-----------+-----------------+
| bsmap    | ZS        | +/- (first slot)|
+----------+-----------+-----------------+
