from os.path import join

# Command: identify
DIR_PREDICT = 'predict'
DIR_IDENTIFY_INTERMEDIATE = join(DIR_PREDICT, 'intermediate_results')
DIR_PREDICT_GENES = join(DIR_IDENTIFY_INTERMEDIATE, 'call_genes')
DIR_IDENTIFY_FASTA = join(DIR_IDENTIFY_INTERMEDIATE, 'single_copy_fasta')
PATH_TLN_TABLE_SUMMARY = join(DIR_PREDICT, '{prefix}.translation_table_summary.tsv')
PATH_FEATURES_SUMMARY = join(DIR_PREDICT, '{prefix}.feature_summary.tsv')
PATH_FAILS = join(DIR_PREDICT,'{prefix}.failed_genomes.tsv')

# Command: identify -> marker genes
GENOME_FILE_SUFFIX = "_genomic.fna"
PROTEIN_FILE_SUFFIX = "_protein.faa"
NT_GENE_FILE_SUFFIX = "_protein.fna"
GFF_FILE_SUFFIX = "_protein.gff"
TRANSLATION_TABLE_SUFFIX = "_translation_table.tsv"
CHECKSUM_SUFFIX = ".sha256"


# General files
PATH_WARNINGS = '{prefix}.warnings.log'