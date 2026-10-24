.. _files/feature_summary.tsv:

feature_summary.tsv
===================

Tab-separated file containing the eight features used by the classifiers for each genome.

* ``detect_table`` writes it to ``predict/[prefix].feature_summary.tsv`` (default prefix: ``gtranslate``),
  with a symbolic link of the same name in the output directory.
* ``build_features`` writes it to ``gtranslate_training.feature_summary.tsv`` in the output directory. This
  is the file to provide to ``fit_models``.

Genomes skipped because Prodigal could not call genes are not included (see
:ref:`failed_genomes.tsv <files/failed_genomes.tsv>`). Values are rounded to 5 decimal places.

The columns in this file are as follows:

* ``user_genome`` (``Genome ID`` for ``build_features``): unique identifier of the genome.
* ``predicted_tln_table`` (``detect_table`` only): the predicted translation table (4, 11 or 25).
* ``Coding_density_4``: gene coding density (%) when predicting genes with Prodigal using translation table 4,
  or equivalently 25, which reassign the UGA stop codon to tryptophan or glycine, respectively.
* ``Coding_density_11``: gene coding density (%) when predicting genes with Prodigal (Hyatt et al., 2010) using
  translation table 11.
* ``Density_Diff``: the difference in coding density between translation table 4/25 and translation table 11,
  i.e. ``Coding_density_4`` - ``Coding_density_11``.
* ``GC``: percentage of guanine (G) and cytosine (C) nucleotides in the genome.
* ``Trp_ratio``: natural log of (UGA + 1) / (UGG + 1), where UGA and UGG are the number of in-frame codons in genes
  predicted under translation table 4. Clipped between -6 and 5 to remove extreme outliers.
* ``Trp_magnitude``: natural log of (UGA + UGG + 1) in genes predicted under translation table 4.
* ``Gly_ratio``: natural log of (UGA + 1) / (GGN + 1), where GGN is the number of glycine codons in genes
  predicted under translation table 4. Clipped between -10 and 0 to remove extreme outliers.
* ``UGG_density``: ratio of UGG tryptophan codons to GGN glycine codons in genes predicted under translation
  table 4.


Produced by
-----------

* :ref:`commands/build_features`
* :ref:`commands/detect_table`
