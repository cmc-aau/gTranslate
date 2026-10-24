.. _about/how_it_works:

How gTranslate works
====================

gTranslate predicts whether a prokaryotic genome uses translation table 11 (the standard bacterial code), 4 (UGA
reassigned from a stop codon to tryptophan) or 25 (UGA reassigned to glycine). It does this in three steps:
calling genes under different translation tables, summarising the gene calls as eight features, and predicting the
translation table with an ensemble of five machine learning classifiers.

.. raw:: html

   <object data="../_static/gtranslate_workflow.svg" type="image/svg+xml"
           style="width:100%;aspect-ratio:960/632" aria-label="How gTranslate predicts a translation table">
     <img src="../_static/gtranslate_workflow.svg" alt="How gTranslate predicts a translation table">
   </object>


1. Calling genes under translation tables 11 and 4
---------------------------------------------------

Genes are called twice with `Prodigal <https://github.com/hyattpd/Prodigal>`_ (Hyatt et al., 2010), once using
translation table 11 and once using translation table 4. Table 4 is used to represent both UGA reassignments, as
tables 4 and 25 only differ in the amino acid encoded by UGA and therefore result in the same gene calls.

If a genome reassigns UGA, genes called under table 11 are truncated at each in-frame UGA codon, which lowers the
coding density. Under table 4, these genes are read through and the coding density is higher. For genomes using
table 11, both tables give similar coding densities.

Prodigal is run with restricted gene prediction across ambiguous bases (``-m``). Genomes under 100 kb are processed
with Prodigal's metagenome mode (``-p meta``). Larger genomes use its single genome mode, falling back to the
metagenome mode if Prodigal fails, in which case a warning is reported.


2. Building the feature vector
------------------------------

The gene calls are summarised as eight features:

* the coding density under translation tables 4 and 11, and the difference between them;
* the GC content of the genome;
* the ratio of in-frame UGA to UGG (tryptophan) codons and the total number of UGA and UGG codons;
* the ratio of in-frame UGA to GGN (glycine) codons, and the ratio of UGG to GGN codons.

The codon features are calculated from the genes called under table 4. As UGA is read as tryptophan in table 4 and
glycine in table 25, the codon usage around UGA helps distinguish between these two tables. The exact definitions are
given in :ref:`files/feature_summary.tsv`.


3. Predicting the translation table
-----------------------------------

The features are provided to five classifiers: AdaBoost, decision tree, k-nearest neighbours, multi-layer perceptron,
and XGBoost (see :ref:`about/models`). The translation table is predicted by majority vote. Each classifier's vote and
a confidence score are reported in :ref:`files/translation_table_summary.tsv`, along with how ties are resolved.

If table 25 is predicted, genes are called again with Prodigal using table 25, so that the genes kept with
``--keep_called_genes`` are translated using the predicted table.


Limitations
-----------

* gTranslate only predicts translation tables 11, 4 and 25, i.e. whether the UGA stop codon is reassigned to
  tryptophan or glycine. Other reassignments, such as UAG to pyrrolysine, are not predicted.
* gTranslate is designed for prokaryotic genomes and should not be applied to eukaryotic or organelle genomes.
* Predictions are most reliable for genomes similar to those used for training. Predictions for genomes from
  taxa without close relatives in the training data, and for highly incomplete or fragmented assemblies, should be
  treated with more caution (see :ref:`about/validation`).
