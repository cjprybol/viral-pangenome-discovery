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

# how to pass parameters via command line
# snakemake --config muscle-params="-msf"
# {config[muscle-params]}"  
# [--config [KEY=VALUE [KEY=VALUE ...]]]
# [--configfile FILE [FILE ...]]
    
    
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
        "data/taxon_10239.genbank/"
    output:
        "data/taxon_10239.genbank/rehydrated.done"
    resources:
        tmpdir="data/"
    shell:
        """
        datasets rehydrate --directory data/taxon_10239.genbank
        touch data/taxon_10239.genbank/rehydrated.done
        """

# snakemake --snakefile Snakefile.py --cores 1 merge_reference_fastas --dry-run
rule merge_reference_fastas:
    input:
        "data/taxon_10239.genbank"
    output:
        "data/taxon_10239.genbank/joint.fasta"
    shell:
        """
        papermill merge-fastas.ipynb
        """

################################################################################################
# DOWNLOAD DATA
################################################################################################

# download IMG/VR
# https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=IMG_VR
# need globus to do automatically
# don't do that, just download, upload to google drive, and then pull from that as a step


# converting jsonl to tables
# dataformat tsv genome --inputfile human/ncbi_dataset/data/assembly_data_report.jsonl

# snakemake --snakefile Snakefile.py --cores 1 sra_prefetch_mt_pleasant
rule sra_prefetch_mt_pleasant:
    input:
        "metadata/usa_mt-pleasant-research-farm_cornell-university_new-york/SraAccList.txt"
    output:
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/prefetch.done"
    shell:
        """
        prefetch \
            --max-size u \
            --progress \
            --output-directory data/usa_mt-pleasant-research-farm_cornell-university_new-york/ \
            --option-file metadata/usa_mt-pleasant-research-farm_cornell-university_new-york/SraAccList.txt
        touch {output}
        """

# snakemake --snakefile Snakefile.py --cores 1 sra_prefetch_exposome
rule sra_prefetch_exposome:
    input:
        "metadata/exposome/SraAccList-10.txt"
        # "metadata/exposome/SraAccList.txt"
    # output:
    #     "data/exposome/prefetch.done"
    shell:
        """
        prefetch \
            --max-size u \
            --progress \
            --output-directory data/exposome/ \
            --option-file {input}
        rm -r SRR*
        rm GL000*
        rm CM0006*
        touch {output}
        """

# work/viral-pangenome-discovery/metadata//SraAccList-10.txt
# snakemake --snakefile Snakefile.py --cores 1 sra_fasterq_exposome
rule sra_fasterq_exposome_459:
    shell:
        """
        fasterq-dump --split-3 data/exposome/SRR6399459.sralite -O data/exposome/SRR6399459.sralite
        mv data/exposome/SRR6399459/SRR6399459.sralite_1.fastq data/exposome/SRR6399459/SRR6399459_1.fastq
        mv data/exposome/SRR6399459/SRR6399459.sralite_2.fastq data/exposome/SRR6399459/SRR6399459_2.fastq
        gzip data/exposome/SRR6399459/*.fastq
        rm data/exposome/SRR6399459.sralite*
        """


# snakemake --snakefile Snakefile.py --cores 1 sra_prefetch_uc_boulder_wastewater
rule sra_prefetch_uc_boulder_wastewater:
    input:
        "metadata/uc_boulder_wastewater/SraAccList-10.txt"
        # "metadata/uc_boulder_wastewater/SraAccList.txt"
    output:
        "data/uc_boulder_wastewater/prefetch.done"
    shell:
        """
        prefetch \
            --max-size u \
            --progress \
            --output-directory data/uc_boulder_wastewater/ \
            --option-file {input}
        touch {output}
        """

# this did not work as expected
# wrote all of the fastqs to the same output directory
# uc_boulder_data_path='data/uc_boulder_wastewater/'
# # snakemake --snakefile Snakefile.py --cores 1 sra_fasterq_uc_boulder_wastewater
# rule sra_fasterq_uc_boulder_wastewater:
#     input:
#         expand("{srr}", srr=[uc_boulder_data_path + '{}'.format(x) 
#                                  for x in os.listdir(uc_boulder_data_path)
#                                      if (os.path.isdir(uc_boulder_data_path + '{}'.format(x))
#                                          and x.startswith("SRR"))]
#               )
#     shell:
#         """
#         fasterq-dump {input} -O {input} --split-3
#         gzip {input}/*.fastq
#         """
        
# rule all:
#     input
        
# snakemake --snakefile Snakefile.py --cores 1 sra_fasterq_uc_boulder_wastewater
# rule sra_fasterq_uc_boulder_wastewater:
#     input:
#         'data/uc_boulder_wastewater/{SRR}/fasterq.done'
#     output:
#         'data/uc_boulder_wastewater/{SRR}/fasterq
        
    
# human microbiome data
# snakemake --snakefile Snakefile.py --cores 1 download_hmp_data
rule download_hmp_data:
    shell:
        """
        echo "implement me"
        """

################################################################################################
# SUBSAMPLE READS FOR DEVELOPMENT TESTING
################################################################################################

# snakemake --snakefile snakemake/0.data-acquisition.snakemake.py --cores all subsample_SRR6399459_10k
rule subsample_SRR6399459_10k:
    shell:
        """
        seqtk sample -s100 data/SRA/SRR6399459/SRR6399459_1.fastq.gz 0.0001 | gzip -c > data/SRA/SRR6399459/SRR6399459_1.10k.fastq.gz
        seqtk sample -s100 data/SRA/SRR6399459/SRR6399459_2.fastq.gz 0.0001 | gzip -c > data/SRA/SRR6399459/SRR6399459_2.10k.fastq.gz
        """

# snakemake --snakefile snakemake/0.data-acquisition.snakemake.py --cores all subsample_SRR6399459_100k
rule subsample_SRR6399459_100k:
    shell:
        """
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_1.fastq.gz 0.001 | gzip -c > data/exposome/SRR6399459/SRR6399459_1.100k.fastq.gz
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_2.fastq.gz 0.001 | gzip -c > data/exposome/SRR6399459/SRR6399459_2.100k.fastq.gz
        """
        
# snakemake --snakefile snakemake/0.data-acquisition.snakemake.py --cores all subsample_SRR6399459_1M
rule subsample_SRR6399459_1M:
    shell:
        """
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_1.fastq.gz 0.01 | gzip -c > data/exposome/SRR6399459/SRR6399459_1.1M.fastq.gz
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_2.fastq.gz 0.01 | gzip -c > data/exposome/SRR6399459/SRR6399459_2.1M.fastq.gz
        """
        
# snakemake --snakefile snakemake/0.data-acquisition.snakemake.py --cores all subsample_SRR6399459_10M
rule subsample_SRR6399459_10M:
    shell:
        """
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_1.fastq.gz 0.1 | gzip -c > data/exposome/SRR6399459/SRR6399459_1.10M.fastq.gz
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_2.fastq.gz 0.1 | gzip -c > data/exposome/SRR6399459/SRR6399459_2.10M.fastq.gz
        """