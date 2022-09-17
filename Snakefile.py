# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# snakemake --lint

import os


# install centrifuge
# git clone https://github.com/infphilo/centrifuge
# cd centrifuge
# make
# # add this to my path
# make install prefix=/home/jovyan/.local

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

# snakemake --snakefile Snakefile.py --cores 1 sra_fasterq_mt_pleasant
rule sra_fasterq:
    input:
        # [ for dataset in DATASETS]
        # ['data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x) for x in os.listdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york') if (os.path.isdir('data/usa_mt-pleasant-research-farm_cornell-university_new-york/{}'.format(x)) and x.startswith("SRR"))]
        # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469"
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476619"
        # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/prefetch.done"
    output:
        # "data/usa_mt-pleasant-research-farm_cornell-university_new-york/fasterq.done"
    shell:
        """
        fasterq-dump --split-3 {input} -O {input}
        gzip {input}/*.fastq
        """
        
# snakemake --snakefile Snakefile.py --cores 1 sra_prefetch_exposome
rule sra_prefetch_exposome:
    input:
        "metadata/exposome/SraAccList-10.txt"
        # "metadata/exposome/SraAccList.txt"
    output:
        "data/exposome/prefetch.done"
    shell:
        """
        prefetch \
            --max-size u \
            --progress \
            --output-directory data/exposome/ \
            --option-file {input}
        touch {output}
        """

# work/viral-pangenome-discovery/metadata//SraAccList-10.txt

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

# human microbiome data
# snakemake --snakefile Snakefile.py --cores 1 download_hmp_data
rule download_hmp_data:
    shell:
        """
        echo "implement me"
        """

################################################################################################
# GET INDICES
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

# https://benlangmead.github.io/aws-indexes/centrifuge
# NCBI: nucleotide non-redundant sequences 	March, 2018
#
# https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
# snakemake --snakefile Snakefile.py --cores 1 download_centrifuge_nt
# done!
rule download_centrifuge_nt:
    output:
        directory("data/centrifuge-databases/nt_2018_3_3")
    shell:
        """
        mkdir -p {output}
        wget --no-clobber --directory-prefix data/centrifuge-databases/ https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
        tar -xvzf data/centrifuge-databases/nt_2018_3_3.tar.gz --directory {output}
        """
        
################################################################################################
# CLASSIFY READS
################################################################################################

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

# Loading database information... done.
# Processed 1790000 sequences (540580000 bp) ...
# 106359342 sequences (32120.52 Mbp) processed in 2731.569s (2336.2 Kseq/m, 705.54 Mbp/m).
#   9757084 sequences classified (9.17%)
#   96602258 sequences unclassified (90.83%)
# snakemake --snakefile Snakefile.py --cores 1 classify_kraken_pluspfp16_db_mt_pleasant
rule classify_kraken_pluspfp16_db_mt_pleasant:
    input: 
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz", 
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz"
    output:
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_16gb_20220607/kraken-output.txt"
    shell:
        """
        kraken2 \
            --report-zero-counts \
            --use-names \
            --threads `nproc` \
            --db "data/kraken-databases/k2_pluspfp_16gb_20220607" \
            --output {output} \
            --report "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_16gb_20220607/kraken-report.txt" \
            --gzip-compressed \
            --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_16gb_20220607/SRR6476469#.classified.fastq" \
            --paired {input}
        """
        
# ????????????
# TODO
# snakemake --snakefile Snakefile.py --cores 1 classify_kraken_pluspfp_db_mt_pleasant
rule classify_kraken_pluspfp_db_mt_pleasant:
    input: 
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz", 
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz"
    output:
        "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_20220607/kraken-output.txt"
    shell:
        """
        kraken2 \
            --memory-mapping \
            --report-zero-counts \
            --use-names \
            --threads `nproc` \
            --db "data/kraken-databases/k2_pluspfp_20220607" \
            --output {output} \
            --report "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_20220607/kraken-report.txt" \
            --gzip-compressed \
            --classified-out "data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/k2_pluspfp_20220607/SRR6476469#.classified.fastq" \
            --paired {input}
        """
        
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
    
# rule trimmomatic:

# rule index_fasta:
#     input:
#     output:
#     samtools faidx:


# Total time for call to driver() for forward index: 00:48:53
# hisat2-build -p `nproc` data/taxon_10239.genbank/joint.fna data/taxon_10239.genbank/joint.fna
# hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
# --al <path>, --al-gz <path>, --al-bz2 <path>

# --al-conc-gz /path/to/file_prefix%.aligned.fastq.gz
# --no-unal
# -p/--threads

# Total time for backward call to driver() for mirror index: 00:18:55
# bowtie2-build --threads `nproc` data/taxon_10239.genbank/joint.fna data/taxon_10239.genbank/joint.fna


# bwa index data/taxon_10239.genbank/joint.fna


# need to pass through samtools to filter out unaligned & secondary
# bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
# now try mapping with each of them

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

# # https://benlangmead.github.io/aws-indexes/centrifuge
# # NCBI: nucleotide non-redundant sequences 	March, 2018
# #
# https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz


########################################################################################
# assemble into contigs
########################################################################################


########################################################################################
# classify assembled contigs
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
        cd blastdb
        time update_blastdb.pl --decompress nt
        """

# update_blastdb.pl --showall pretty
# refseq_protein for protein
# nt for DNA (or env_nt for even more metagenomic projects)
# time update_blastdb.pl --decompress nt
# https://github.com/ncbi/blast_plus_docs#blast-databases
# nt