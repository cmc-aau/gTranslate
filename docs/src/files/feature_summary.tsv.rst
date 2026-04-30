.. _files/feature_summary.tsv:

feature_summary.tsv
====================

The tab-delimited file containing the summary of features for each genome in the training dataset.

The columns in this file are as follows:
* user_genome ( or genome ID for ``build_features``): Unique identifier of the query genome.
* best_tln_table ( missing for ``build_features``): The final predicted Genetic Translation Table (GTT). This is the "consensus" choice recommended by the tool (e.g., 11 for standard, 4 for UGA=Trp, or 25 for UGA=Gly).
* Coding_density_4: gene coding density when predicting genes with Prodigal using translation table 4 or equivalently 25 which reassigned the UGA stop coding to either tryptophan or glycine.
* Coding_density_11: gene coding density when predicting genes with Prodigal (Hyatt et al., 2010) using translation table 11.
* Density_Diff: the difference in coding density when using translation table 4/25 or translation table 11, i.e. CD4 – CD11.
* GC: percentage of guanine (G) and cytosine (C) nucleotides in a genome.
* Trp_ratio: log-transformed ratio of UGA to UGG codon counts when predicting genes with Prodigal under translation table 4. The log ratio is clipped between -6 and 5 to remove extreme outliers.
* Trp_magnitude: log-transformed count of all UGA and UGG tryptophan codons when predicting genes with Prodigal under translation table 4.
* Gly_ratio: log-transformed ratio of UGA codon counts to the total glycine codon counts (i.e. codons GGn) when predicting genes with Prodigal under translation table 4. The log ratio is clipped between -10 and 0 to remove extreme outliers.
* UGG_density: ratio of UGG tryptophan codons to glycine codons (GGn) when predicting genes with Prodigal under translation table 4.


Produced by
-----------

* :ref:`commands/build_features`
* :ref:`commands/detect_table`