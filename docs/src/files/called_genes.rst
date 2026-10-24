.. _files/called_genes:

Called genes
============

Genes called by Prodigal under the predicted translation table. These are written to
``predict/intermediate_results/call_genes/[genome_id]/`` and deleted at the end of the run unless
``--keep_called_genes`` is specified.

* ``[genome_id]_protein.faa``: protein sequences of the called genes in FASTA format.
* ``[genome_id]_protein.fna``: nucleotide sequences of the called genes in FASTA format.
* ``[genome_id]_protein.gff``: coordinates of the called genes in GFF format.
* ``[genome_id]_translation_table.tsv``: the prediction, confidence, warnings and features for the genome, as
  reported in :ref:`translation_table_summary.tsv <files/translation_table_summary.tsv>`.
* ``*.sha256``: checksums of the above files, used to skip genomes that have already been processed when
  ``detect_table`` is re-run with the same output directory.


Produced by
-----------

* :ref:`commands/detect_table` (with ``--keep_called_genes``)
