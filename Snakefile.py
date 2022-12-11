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

# input: "path/to/inputfile", "path/to/other/inputfile"
# output: "path/to/outputfile", somename = "path/to/another/outputfile"
# run:
#     for f in input:
#         ...
#         with open(output[0], "w") as out:
#             out.write(...)
#     with open(output.somename, "w") as out:
#         out.write(...)
        
# # snakemake --snakefile Snakefile.py --cores 1 sra_fasterq_mt_pleasant
# rule sra_fasterq_mt_pleasant:
#     input:
#         # [ for dataset in DATASETS]
#         # ['data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x) for x in os.listdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york') if (os.path.isdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x)) and x.startswith("SRR"))]
#         # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469"
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476619"
#         # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/prefetch.done"
#     output:
#         # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/fasterq.done"
#     shell:
#         """
#         fasterq-dump --split-3 {input} -O {input}
#         gzip {input}/*.fastq
#         """
        
# input: "path/to/inputfile", "path/to/other/inputfile"
# output: "path/to/outputfile", somename = "path/to/another/outputfile"
# run:
#     for f in input:
#         ...
#         with open(output[0], "w") as out:
#             out.write(...)
#     with open(output.somename, "w") as out:
#         out.write(...)
    
    
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
# SUBSAMPLE READS
################################################################################################

# snakemake --snakefile Snakefile.py --cores all subsample_SRR6399459_100k
rule subsample_SRR6399459_100k:
    shell:
        """
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_1.fastq.gz 0.001 | gzip -c > data/exposome/SRR6399459/SRR6399459_1.100k.fastq.gz
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_2.fastq.gz 0.001 | gzip -c > data/exposome/SRR6399459/SRR6399459_2.100k.fastq.gz
        """
        
# snakemake --snakefile Snakefile.py --cores all subsample_SRR6399459_1M
rule subsample_SRR6399459_1M:
    shell:
        """
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_1.fastq.gz 0.01 | gzip -c > data/exposome/SRR6399459/SRR6399459_1.1M.fastq.gz
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_2.fastq.gz 0.01 | gzip -c > data/exposome/SRR6399459/SRR6399459_2.1M.fastq.gz
        """
        
# snakemake --snakefile Snakefile.py --cores all subsample_SRR6399459_10M
rule subsample_SRR6399459_10M:
    shell:
        """
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_1.fastq.gz 0.1 | gzip -c > data/exposome/SRR6399459/SRR6399459_1.10M.fastq.gz
        seqtk sample -s100 data/exposome/SRR6399459/SRR6399459_2.fastq.gz 0.1 | gzip -c > data/exposome/SRR6399459/SRR6399459_2.10M.fastq.gz
        """

################################################################################################
# QC READS (fastqc and trimmomatic)
################################################################################################

# rule trimmomatic:
        
################################################################################################
# CLASSIFY READS WITH KRAKEN
################################################################################################

# snakemake --snakefile Snakefile.py --cores 1 download_kraken_viruses_db
# rule download_kraken_viruses_db:
#     output:
#         directory("data/kraken-databases/k2_viral_20220607")
#     shell:
#         """
#         wget --no-clobber --directory-prefix data/kraken-databases/ https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220607.tar.gz
#         mkdir -p {output}
#         tar -xvzf data/kraken-databases/k2_viral_20220607.tar.gz
#         """
# archaea, bacteria, viral, plasmid, human1, UniVec_Core
# Standard plus protozoa & fungi
# Standard plus protozoa, fungi & plant
# capped at 8
# snakemake --snakefile Snakefile.py --cores 1 download_kraken_pluspfp8_db
# rule download_kraken_pluspfp8_db:
#     output:
#         directory("data/kraken-databases/k2_pluspfp_08gb_20220607")
#     shell:
#         """
#         wget --no-clobber --directory-prefix data/kraken-databases/ https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20220607.tar.gz
#         mkdir -p {output}
#         tar -xvzf data/kraken-databases/k2_pluspfp_08gb_20220607.tar.gz --directory {output}
#         """

# done!
# capped at 16
# snakemake --snakefile Snakefile.py --cores 1 download_kraken_pluspfp16_db
rule download_kraken_pluspfp16_db:
    output:
        directory("data/kraken-databases/k2_pluspfp_16gb_20220607")
    shell:
        """
        wget --no-clobber --directory-prefix data/kraken-databases/ https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20220607.tar.gz
        mkdir -p {output}
        tar -xvzf data/kraken-databases/k2_pluspfp_16gb_20220607.tar.gz --directory {output}
        """
# TODO
# snakemake --snakefile Snakefile.py --cores 1 download_kraken_pluspfp_db
rule download_kraken_pluspfp_db:
    output:
        directory("data/kraken-databases/k2_pluspfp_20220607")
    shell:
        """
        wget --no-clobber --directory-prefix data/kraken-databases/ https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20220607.tar.gz
        mkdir -p {output}
        tar -xvzf data/kraken-databases/k2_pluspfp_20220607.tar.gz --directory {output}
        """

# Loading database information... done.
# 106359342 sequences (32120.52 Mbp) processed in 802.851s (7948.6 Kseq/m, 2400.48 Mbp/m).
#   150693 sequences classified (0.14%)
#   106208649 sequences unclassified (99.86%)
# # snakemake --snakefile Snakefile.py --cores 1 download_kraken_viruses_db
# rule classify_kraken_viruses_db_mt_pleasant:
#     input:
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz"
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz"
#         # ['data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x) for x in os.listdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york') if (os.path.isdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x)) and x.startswith("SRR"))]
#         # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469"
#         # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/prefetch.done"
#     db:
#         "data/kraken-databases/k2_viral_20220607"
#     output:
#         output="data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-output.txt"
#         report="data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-report.txt"
#     shell:
#         """
#         kraken2 \
#             --report-zero-counts \
#             --use-names \
#             --threads `nproc` \
#             --db {db} \
#             --output {output.output} \
#             --report {output.report} \
#             --gzip-compressed \
#             --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/SRR6476469#.classified.fastq" \
#             --paired {input}
#         """


# ['data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x) for x in os.listdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york') if (os.path.isdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x)) and x.startswith("SRR"))]
# "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469"
# "data/usa_mt-pleasant-research-farm_cornell-university_new-york/prefetch.done"


# classified way more sequences using the 16gb than the 8
# # Loading database information... done.
# # Processed 5410000 sequences (1633820000 bp) ...
# # 106359342 sequences (32120.52 Mbp) processed in 2363.352s (2700.2 Kseq/m, 815.47 Mbp/m).
# #   5348645 sequences classified (5.03%)
# #   101010697 sequences unclassified (94.97%)

# # snakemake --snakefile Snakefile.py --cores 1 classify_kraken_pluspfp8_db_mt_pleasant
# rule classify_kraken_pluspfp8_db_mt_pleasant:
#     input: 
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz", 
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz"
#     output:
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_08gb_20220607/kraken-output.txt"
#     shell:
#         """
#         kraken2 \
#             --report-zero-counts \
#             --use-names \
#             --threads `nproc` \
#             --db "data/kraken-databases/k2_pluspfp_08gb_20220607" \
#             --output {output} \
#             --report "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_08gb_20220607/kraken-report.txt" \
#             --gzip-compressed \
#             --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_08gb_20220607/SRR6476469#.classified.fastq" \
#             --paired {input}
#         """

# # Loading database information... done.
# # Processed 1790000 sequences (540580000 bp) ...
# # 106359342 sequences (32120.52 Mbp) processed in 2731.569s (2336.2 Kseq/m, 705.54 Mbp/m).
# #   9757084 sequences classified (9.17%)
# #   96602258 sequences unclassified (90.83%)
# # snakemake --snakefile Snakefile.py --cores 1 classify_kraken_pluspfp16_db_mt_pleasant
# rule classify_kraken_pluspfp16_db_mt_pleasant:
#     input: 
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz", 
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz"
#     output:
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_16gb_20220607/kraken-output.txt"
#     shell:
#         """
#         kraken2 \
#             --report-zero-counts \
#             --use-names \
#             --threads `nproc` \
#             --db "data/kraken-databases/k2_pluspfp_16gb_20220607" \
#             --output {output} \
#             --report "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_16gb_20220607/kraken-report.txt" \
#             --gzip-compressed \
#             --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_16gb_20220607/SRR6476469#.classified.fastq" \
#             --paired {input}
#         """
        
# # ????????????
# # TODO
# # snakemake --snakefile Snakefile.py --cores 1 classify_kraken_pluspfp_db_mt_pleasant
# rule classify_kraken_pluspfp_db_mt_pleasant:
#     input: 
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz", 
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz"
#     output:
#         "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_20220607/kraken-output.txt"
#     shell:
#         """
#         kraken2 \
#             --memory-mapping \
#             --report-zero-counts \
#             --use-names \
#             --threads `nproc` \
#             --db "data/kraken-databases/k2_pluspfp_20220607" \
#             --output {output} \
#             --report "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_20220607/kraken-report.txt" \
#             --gzip-compressed \
#             --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_20220607/SRR6476469#.classified.fastq" \
#             --paired {input}
#         """

# Running!!
# snakemake --snakefile Snakefile.py --cores all classify_kraken_pluspfp_db_SRR6399459_100k
rule classify_kraken_pluspfp16_SRR6399459_100k:
    input: 
        "data/exposome/SRR6399459/SRR6399459_1.fastq.gz", 
        "data/exposome/SRR6399459/SRR6399459_2.fastq.gz"
    output:
        "data/exposome/SRR6399459/k2_pluspfp_16gb_20220607/kraken-output.txt"
    shell:
        """
        kraken2 \
            --report-zero-counts \
            --use-names \
            --threads `nproc` \
            --db "data/kraken-databases/k2_pluspfp_16gb_20220607" \
            --output {output} \
            --report "data/exposome/SRR6399459//k2_pluspfp_16gb_20220607/kraken-report.txt" \
            --gzip-compressed \
            --classified-out "data/exposome/SRR6399459/k2_pluspfp_16gb_20220607/SRR6476469#.classified.fastq" \
            --paired {input}
        """
        
# # snakemake --snakefile Snakefile.py --cores all classify_kraken_pluspfp_db_SRR6399459_1M
# rule classify_kraken_pluspfp16_SRR6399459_1M:
#     input: 
#         "data/exposome/SRR6399459/SRR6399459_1.fastq.gz", 
#         "data/exposome/SRR6399459/SRR6399459_2.fastq.gz"
#     output:
#         "data/exposome/SRR6399459/k2_pluspfp_16gb_20220607/kraken-output.txt"
#     shell:
#         """
#         kraken2 \
#             --report-zero-counts \
#             --use-names \
#             --threads `nproc` \
#             --db "data/kraken-databases/k2_pluspfp_16gb_20220607" \
#             --output {output} \
#             --report "data/exposome/SRR6399459//k2_pluspfp_16gb_20220607/kraken-report.txt" \
#             --gzip-compressed \
#             --classified-out "data/exposome/SRR6399459/k2_pluspfp_16gb_20220607/SRR6476469#.classified.fastq" \
#             --paired {input}
#         """
        
# # snakemake --snakefile Snakefile.py --cores all classify_kraken_pluspfp_db_SRR6399459_10M
# rule classify_kraken_pluspfp16_SRR6399459_10M:
#     input: 
#         "data/exposome/SRR6399459/SRR6399459_1.fastq.gz", 
#         "data/exposome/SRR6399459/SRR6399459_2.fastq.gz"
#     output:
#         "data/exposome/SRR6399459/k2_pluspfp_16gb_20220607/kraken-output.txt"
#     shell:
#         """
#         kraken2 \
#             --report-zero-counts \
#             --use-names \
#             --threads `nproc` \
#             --db "data/kraken-databases/k2_pluspfp_16gb_20220607" \
#             --output {output} \
#             --report "data/exposome/SRR6399459//k2_pluspfp_16gb_20220607/kraken-report.txt" \
#             --gzip-compressed \
#             --classified-out "data/exposome/SRR6399459/k2_pluspfp_16gb_20220607/SRR6476469#.classified.fastq" \
#             --paired {input}
#         """
        
# snakemake --snakefile Snakefile.py --cores all classify_kraken_pluspfp_db_SRR6399459_100k
rule classify_kraken_pluspfp_db_SRR6399459_100k:
    input: 
        "data/exposome/SRR6399459/SRR6399459_1.100k.fastq.gz", 
        "data/exposome/SRR6399459/SRR6399459_2.100k.fastq.gz"
    output:
        "data/exposome/SRR6399459/k2_pluspfp_20220607/SRR6399459_1.100k.fastq.kraken-output.txt"
    shell:
        """
        kraken2 \
            --memory-mapping \
            --report-zero-counts \
            --use-names \
            --threads `nproc` \
            --db "data/kraken-databases/k2_pluspfp_20220607" \
            --output {output} \
            --report "data/exposome/SRR6399459/k2_pluspfp_20220607/k2_pluspfp_20220607/SRR6399459_1.100k.fastq.kraken-report.txt" \
            --gzip-compressed \
            --classified-out "data/exposome/SRR6399459/k2_pluspfp_20220607/k2_pluspfp_20220607/SRR6476469#.100k.classified.fastq" \
            --paired {input}
        """
    
# # snakemake --snakefile Snakefile.py --cores all classify_kraken_pluspfp_db_SRR6399459
# rule classify_kraken_pluspfp_db_SRR6399459:
#     input: 
#         "data/exposome/SRR6399459/SRR6399459_1.fastq.gz", 
#         "data/exposome/SRR6399459/SRR6399459_2.fastq.gz"
#     output:
#         "data/exposome/SRR6399459/k2_pluspfp_20220607/kraken-output.txt"
#     shell:
#         """
#         kraken2 \
#             --memory-mapping \
#             --report-zero-counts \
#             --use-names \
#             --threads `nproc` \
#             --db "data/kraken-databases/k2_pluspfp_20220607" \
#             --output {output} \
#             --report "data/exposome/SRR6399459/k2_pluspfp_20220607/k2_pluspfp_20220607/kraken-report.txt" \
#             --gzip-compressed \
#             --classified-out "data/exposome/SRR6399459/k2_pluspfp_20220607/k2_pluspfp_20220607/SRR6476469#.classified.fastq" \
#             --paired {input}
#         """
        
# kraken2 \
#     --report-zero-counts \
#     --use-names \
#     --threads `nproc` \
#     --db data/kraken-databases/k2_viral_20220607 \
#     --output data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-output.txt \
#     --report data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-report.txt \
#     --gzip-compressed \
#     --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/SRR6476469#.classified.fastq" \
#     --paired data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz


# mkdir -p "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607"
# kraken2 \
#     --report-zero-counts \
#     --use-names \
#     --threads `nproc` \
#     --db data/kraken-databases/k2_viral_20220607 \
#     --output data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-output.txt \
#     --report data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-report.txt \
#     --gzip-compressed \
#     --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/SRR6476469#.classified.fastq" \
#     --paired data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz

# wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220607.tar.gz
# tar -xvzf k2_viral_20220607.tar.gz

# didn't create a new directory and didn't actually write out any of the outputs or anything new??

# output file doesn't seem to be helpful at all
# also very very large
# consider throwing away?

# mkdir -p "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607"
# kraken2 \
#     --report-zero-counts \
#     --use-names \
#     --threads `nproc` \
#     --db data/kraken-databases/k2_viral_20220607 \
#     --output data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-output.txt \
#     --report data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/kraken-report.txt \
#     --gzip-compressed \
#     --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_viral_20220607/SRR6476469#.classified.fastq" \
#     --paired data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz

# # archaea, bacteria, viral, plasmid, human1, UniVec_Core
# # Standard plus protozoa & fungi
# # Standard plus protozoa, fungi & plant
# # capped at 8
# https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20220607.tar.gz
# # capped at 16
# https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20220607.tar.gz

########################################################################################
# classify reads with centriguge
########################################################################################

# # https://benlangmead.github.io/aws-indexes/centrifuge
# # NCBI: nucleotide non-redundant sequences 	March, 2018
# #
# https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz

# couldn't get this to run
# https://benlangmead.github.io/aws-indexes/centrifuge
# NCBI: nucleotide non-redundant sequences 	March, 2018
#
# https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
# snakemake --snakefile Snakefile.py --cores 1 download_centrifuge_nt
# done!
# rule download_centrifuge_nt:
#     output:
#         directory("data/centrifuge-databases/nt_2018_3_3")
#     shell:
#         """
#         mkdir -p {output}
#         wget --no-clobber --directory-prefix data/centrifuge-databases/ https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
#         tar -xvzf data/centrifuge-databases/nt_2018_3_3.tar.gz --directory {output}
#         """

# centrifuge [options]* -x <cf-idx> {-1 <m1> -2 <m2> | -U <r> | --sample-sheet <s> } [-S <filename>] [--report-file <report>]

#   -p/--threads <int> number of alignment threads to launch (1)
#   <cf-idx>   Index filename prefix (minus trailing .X.cf).
#   <m1>       Files with #1 mates, paired with files in <m2>.
#              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#   <m2>       Files with #2 mates, paired with files in <m1>.
#              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#   <r>        Files with unpaired reads.
#              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#   <s>        A TSV file where each line represents a sample.
#   <filename>      File for classification output (default: stdout)
#   <report>   File for tabular report output (default: centrifuge_report.tsv)

#   <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
#   specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.

# mkdir -p data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge
# centrifuge \
#     --threads `nproc` \
#     -x data/centrifuge-databases/nt_2018_3_3/nt \
#     -1 data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz \
#     -2 data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz \
#     --al-conc-gz data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/SRR6476469.aligned.fastq.gz \
#     -S data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/classification.txt \
#     --report-file data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/centrifuge_report.tsv

# install centrifuge
# git clone https://github.com/infphilo/centrifuge
# cd centrifuge
# make
# # add this to my path
# make install prefix=/home/jovyan/.local
    

########################################################################################
# classify reads with bwa
########################################################################################

# rule index_fasta:
#     input:
#     output:
#     samtools faidx:

# bwa index data/taxon_10239.genbank/joint.fna

# need to pass through samtools to filter out unaligned & secondary
# bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
# now try mapping with each of them

# workspace/viral-pangenome-discovery/data/taxon_10239.genbank/joint.fna

########################################################################################
# classify reads with bowtie2/hisat
########################################################################################

# Total time for call to driver() for forward index: 00:48:53
# hisat2-build -p `nproc` data/taxon_10239.genbank/joint.fna data/taxon_10239.genbank/joint.fna
# hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
# --al <path>, --al-gz <path>, --al-bz2 <path>

# --al-conc-gz /path/to/file_prefix%.aligned.fastq.gz
# --no-unal
# -p/--threads

# Total time for backward call to driver() for mirror index: 00:18:55
# bowtie2-build --threads `nproc` data/taxon_10239.genbank/joint.fna data/taxon_10239.genbank/joint.fna


########################################################################################
# assemble into contigs with megahit
########################################################################################

# snakemake --snakefile Snakefile.py --cores all megahit_assemble_mt_pleasant_SRR6476469
# rule megahit_assemble_mt_pleasant_SRR6476469:
#     shell:
#         """
#         megahit \
#             -1 /home/jovyan/workspace/viral-pangenome-discovery/data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz \
#             -2 /home/jovyan/workspace/viral-pangenome-discovery/data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz \
#             -o /home/jovyan/workspace/viral-pangenome-discovery/data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/megahit
#         """

# 2022-12-11 14:06:13 - MEGAHIT v1.2.9
# 2022-12-11 14:06:13 - Using megahit_core with POPCNT and BMI2 support
# 2022-12-11 14:06:13 - Convert reads to binary library
# 2022-12-11 14:06:13 - b'INFO  sequence/io/sequence_lib.cpp  :   75 - Lib 0 (/home/jovyan/workspace/viral-pangenome-discovery/data/exposome/SRR6399459/SRR6399459_1.100k.fastq.gz,/home/jovyan/workspace/viral-pangenome-discovery/data/exposome/SRR6399459/SRR6399459_2.100k.fastq.gz): pe, 114394 reads, 151 max length'
# 2022-12-11 14:06:13 - b'INFO  utils/utils.h                 :  152 - Real: 0.2311\tuser: 0.1197\tsys: 0.0179\tmaxrss: 14940'
# 2022-12-11 14:06:13 - k-max reset to: 141 
# 2022-12-11 14:06:13 - Start assembly. Number of CPU threads 48 
# 2022-12-11 14:06:13 - k list: 21,29,39,59,79,99,119,141 
# 2022-12-11 14:06:13 - Memory used: 121546616832
# 2022-12-11 14:06:13 - Extract solid (k+1)-mers for k = 21 
# 2022-12-11 14:06:16 - Build graph for k = 21 
# 2022-12-11 14:06:19 - Assemble contigs from SdBG for k = 21
# 2022-12-11 14:06:20 - Local assembly for k = 21
# 2022-12-11 14:06:21 - Extract iterative edges from k = 21 to 29 
# 2022-12-11 14:06:21 - Build graph for k = 29 
# 2022-12-11 14:06:23 - Assemble contigs from SdBG for k = 29
# 2022-12-11 14:06:24 - Local assembly for k = 29
# 2022-12-11 14:06:25 - Extract iterative edges from k = 29 to 39 
# 2022-12-11 14:06:25 - Build graph for k = 39 
# 2022-12-11 14:06:27 - Assemble contigs from SdBG for k = 39
# 2022-12-11 14:06:28 - Local assembly for k = 39
# 2022-12-11 14:06:29 - Extract iterative edges from k = 39 to 59 
# 2022-12-11 14:06:29 - Build graph for k = 59 
# 2022-12-11 14:06:31 - Assemble contigs from SdBG for k = 59
# 2022-12-11 14:06:32 - Local assembly for k = 59
# 2022-12-11 14:06:32 - Extract iterative edges from k = 59 to 79 
# 2022-12-11 14:06:33 - Build graph for k = 79 
# 2022-12-11 14:06:35 - Assemble contigs from SdBG for k = 79
# 2022-12-11 14:06:36 - Local assembly for k = 79
# 2022-12-11 14:06:37 - Extract iterative edges from k = 79 to 99 
# 2022-12-11 14:06:37 - Build graph for k = 99 
# 2022-12-11 14:06:39 - Assemble contigs from SdBG for k = 99
# 2022-12-11 14:06:40 - Local assembly for k = 99
# 2022-12-11 14:06:41 - Extract iterative edges from k = 99 to 119 
# 2022-12-11 14:06:41 - Build graph for k = 119 
# 2022-12-11 14:06:43 - Assemble contigs from SdBG for k = 119
# 2022-12-11 14:06:44 - Local assembly for k = 119
# 2022-12-11 14:06:45 - Extract iterative edges from k = 119 to 141 
# 2022-12-11 14:06:45 - Build graph for k = 141 
# 2022-12-11 14:06:47 - Assemble contigs from SdBG for k = 141
# 2022-12-11 14:06:47 - Merging to output final contigs 
# 2022-12-11 14:06:47 - 256 contigs, total 326097 bp, min 202 bp, max 16474 bp, avg 1273 bp, N50 2698 bp
# 2022-12-11 14:06:47 - ALL DONE. Time elapsed: 34.613196 seconds 

# snakemake --snakefile Snakefile.py --cores all megahit_assemble_SRR6399459_100k
rule megahit_assemble_SRR6399459_100k:
    shell:
        """
        megahit \
            -1 /home/jovyan/workspace/viral-pangenome-discovery/data/exposome/SRR6399459/SRR6399459_1.100k.fastq.gz \
            -2 /home/jovyan/workspace/viral-pangenome-discovery/data/exposome/SRR6399459/SRR6399459_2.100k.fastq.gz \
            -o /home/jovyan/workspace/viral-pangenome-discovery/data/exposome/SRR6399459/megahit_100k
        """

########################################################################################
# assemble into contigs with spades
########################################################################################

            # --metaviral \
            # --metaplasmid \
        
# Usage: spades.py [options] -o <output_dir>

# Basic options:
#   -o <output_dir>             directory to store all the resulting files (required)
  
#                               file with trusted contigs
  
#                               file with untrusted contigs

#   -t <int>, --threads <int>   number of threads. [default: 16]
#   -m <int>, --memory <int>    RAM limit for SPAdes in Gb (terminates if exceeded). [default: 250]


# # I couldn't get this to run :(
# # I tried:
# # 3.15.0-4
# # 3.14.0-1
# # 3.13.0
# # snakemake --snakefile Snakefile.py --cores all spades_assemble_mt_pleasant_SRR6476469
# rule spades_assemble_mt_pleasant_SRR6476469:
#     shell:
#         """
#         spades.py \
#             --meta \
#             -1 /home/jovyan/workspace/viral-pangenome-discovery/data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz \
#             -2 /home/jovyan/workspace/viral-pangenome-discovery/data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz \
#             -o /home/jovyan/workspace/viral-pangenome-discovery/data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/spades
#         """

########################################################################################
# classify assembled contigs with blast
########################################################################################

# dropping because of non-inclusive DB
# https://sourmash.readthedocs.io/en/latest/tutorials-lca.html
# curl -L -o genbank-k31.lca.json.gz https://osf.io/4f8n3/download
# sourmash sketch dna -p scaled=1000,k=31 --name-from-first some-genome.fa.gz
# sourmash lca summarize --db genbank-k31.lca.json.gz \
    # --query some-genome.fa.gz.sig

# TODO
# snakemake --snakefile Snakefile.py --cores 1 download_blast_nt
rule download_blast_nt:
    shell:
        """
        mkdir -p data/blastdb
        cd data/blastdb
        time update_blastdb.pl --decompress nt
        """

# update_blastdb.pl --showall pretty
# refseq_protein for protein
# nt for DNA (or env_nt for even more metagenomic projects)
# time update_blastdb.pl --decompress nt
# https://github.com/ncbi/blast_plus_docs#blast-databases
# nt

########################################################################################
# classify assembled contigs with minimap
########################################################################################

########################################################################################
# classify assembled contigs with MMSeq2
########################################################################################


# mmseqs databases UniProtKB/Swiss-Prot data/mmseq2/swissprot tmp

# mmseqs databases UniRef90 data/mmseq2/UniRef90 tmp

# --file-allocation=none

# mmseqs databases UniRef90 data/mmseq2/UniRef90 tmp
# mmseqs createdb examples/DB.fasta targetDB
# mmseqs createtaxdb targetDB tmp
# mmseqs createindex targetDB tmp
# mmseqs easy-taxonomy examples/QUERY.fasta targetDB alnRes tmp