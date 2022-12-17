# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# document
# https://snakemake.readthedocs.io/en/v6.0.3/executing/cli.html#visualization
# snakemake --cores 1 --forceall --dag --snakefile $SNAKEFILE \
# | dot -Tpdf \
# > $SNAKEFILE.dag.pdf

# lint
# snakemake --lint --snakefile $SNAKEFILE

import os
import datetime

invoked_datetime = datetime.datetime.now()
invoked_timestamp = f"{invoked_datetime.year}{invoked_datetime.month}{invoked_datetime.day}{invoked_datetime.hour}{invoked_datetime.second}"
    
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