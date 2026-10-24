*******
gTranslate
*******

**gTranslate** is a machine learning-based command-line tool for predicting the translation table (TT) used by prokaryotic genomes. By analyzing specific sequence features — such as coding density differences, Trp ratios, and Gly ratios — `gTranslate` can accurately distinguish between the genetic codes associated with reassignment of the UGA stop codon, i.e. translation tables 11 (standard prokaryotic code), 4 (UGA=Trp), and 25 (UGA=Gly).

gTranslate is open source and released under the
`GNU General Public License (Version 3) <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

Notifications about gTranslate releases will be available through the GTDB Bluesky account
`<https://bsky.app/profile/ace-gtdb.bsky.social>`_.

Please post questions and issues related to gTranslate on the Issues section of the GitHub repository.


Running gTranslate
==================

1. Install gTranslate and download the models (:ref:`installing`)

2. Access the help documentation :ref:`commands`, or view the program help menu: ``gtranslate -h``

* Note: Individual help can be accessed via the specific command, e.g.: ``gtranslate detect_table -h``


Citing gTranslate
=================
If you use gTranslate, please cite:

    Chaumeil P-A, Hugenholtz P, Parks DH. 2026. `gTranslate: rapid and accurate translation table prediction for prokaryotic genomes <https://doi.org/10.7717/peerj.21825>`_. *PeerJ* 14:e21825. https://doi.org/10.7717/peerj.21825

We also encourage you to cite the third-party dependencies as described in :ref:`references`.


.. toctree::
   :caption: Getting started
   :maxdepth: 1

   announcements
   installing/index
   faq


.. toctree::
   :caption: Running gTranslate
   :maxdepth: 1

   commands/index
   files/index


.. toctree::
   :caption: About
   :maxdepth: 1

   changelog
   references