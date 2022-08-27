# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# snakemake --lint

import os

# conda activate viral-pangenome-discovery

# to specificy conda env as part of the rule
# conda:
#     "environment.yaml"

# https://snakemake.readthedocs.io/en/v6.0.3/executing/cli.html#visualization
# snakemake document --cores 1
rule document:
    output:
        "dag.pdf"
    shell:
        """
        snakemake --forceall --dag | dot -Tpdf > dag.pdf
        """
        
################################################################################################
# BY TAXON ID GENBANK
# Collecting 50,633 genome accessions [===========================>--------------------]  59% 30000/50633
################################################################################################

# snakemake --snakefile Snakefile.py --cores 1 download_taxon_10239_genbank_dehydrated
rule download_taxon_10239_genbank_dehydrated:
    output:
        "data/taxon_10239.genbank.zip"
    resources:
        tmpdir="data/"
    shell:
        """
        datasets download genome taxon 10239 --dehydrated --assembly-source genbank --filename {output}
        """

# snakemake --snakefile Snakefile.py --cores 1 unzip_taxon_10239_genbank
rule unzip_taxon_10239_genbank:
    input:
        "data/taxon_10239.genbank.zip"
    resources:
        tmpdir="data/"
    output:
        directory("data/taxon_10239.genbank")
    shell:
        """
        unzip -d {output} {input}
        """
        
# snakemake --snakefile Snakefile.py --cores 1 rehydrate_taxon_10239_genbank
rule rehydrate_taxon_10239_genbank:
    input:
        "data/taxon_10239.genbank"
    resources:
        tmpdir="data/"
    shell:
        """
        datasets rehydrate --directory data/taxon_10239.genbank
        """

# ################################################################################################
# # BY TAXON ID
# # Collecting 62,332 genome accessions [===============>--------------------------------]  35% 22000/62332
# ################################################################################################

# # snakemake --snakefile Snakefile.py --cores 1 download_taxon_10239_genbank_dehydrated
# rule download_taxon_10239_dehydrated:
#     output:
#         "data/taxon_10239.zip"
#     resources:
#         tmpdir="data/"
#     shell:
#         """
#         datasets download genome taxon 10239 --dehydrated --filename {output}
#         """

# # snakemake --snakefile Snakefile.py --cores 1 unzip_taxon_10239
# rule unzip_taxon_10239:
#     input:
#         "data/taxon_10239.zip"
#     resources:
#         tmpdir="data/"
#     output:
#         directory("data/taxon_10239")
#     shell:
#         """
#         unzip -d {output} {input}
#         """
        
# # snakemake --snakefile Snakefile.py --cores 1 rehydrate_taxon_10239
# rule rehydrate_taxon_10239:
#     input:
#         "data/taxon_10239"
#     resources:
#         tmpdir="data/"
#     shell:
#         """
#         datasets rehydrate --directory data/taxon_10239
#         """

# download IMG/VR
# https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=IMG_VR
# need globus to do automatically
# don't do that, just download, upload to google drive, and then pull from that as a step


# converting jsonl to tables
# dataformat tsv genome --inputfile human/ncbi_dataset/data/assembly_data_report.jsonl