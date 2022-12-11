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
# snakemake --snakefile snakemake/1.fastq-qc.snakemake.py --cores 1 document
rule document:
    output:
        "dag.pdf"
    shell:
        """
        snakemake --forceall --dag | dot -Tpdf > dag.pdf
        """

# [--config [KEY=VALUE [KEY=VALUE ...]]]
# [--configfile FILE [FILE ...]]
    
sample_id = os.path.basename(config["sample"])
# forward_reads = f'{config[sample]}/{sample_id}_1.10k.fastq.gz'
forward_reads = config["sample"] + f'/{sample_id}_1.10k.fastq.gz'
# reverse_reads = f'{config[sample]}/{sample_id}_1.10k.fastq.gz'
reverse_reads = config["sample"] + f'/{sample_id}_2.10k.fastq.gz'
fastqc_outdir = config["sample"] + '/fastqc'

forward_html_report = fastqc_outdir + "/" + sample_id + "_1.10k_fastqc.html"
reverse_html_report = fastqc_outdir + "/" + sample_id + "_2.10k_fastqc.html"
forward_zip = fastqc_outdir + "/" + sample_id + "_1.10k_fastqc.zip"
reverse_zip = fastqc_outdir + "/" + sample_id + "_2.10k_fastqc.zip"

# snakemake --snakefile snakemake/1.fastq-qc.snakemake.py --cores 1 fastqc --config sample=data/SRA/SRR6399459
rule fastqc:
    input: 
        forward_reads,
        reverse_reads
    output: 
        forward_html_report,
        reverse_html_report,
        forward_zip,
        reverse_zip
    shell:
        """
        mkdir -p {fastqc_outdir}
        fastqc --outdir {fastqc_outdir} {forward_reads} {reverse_reads}
        """

trim_galore_outdir = config["sample"] + '/trim_galore'
trimmed_forward_reads = trim_galore_outdir + f'/{sample_id}_1.10k_val_1.fq.gz'
trimmed_reverse_reads = trim_galore_outdir + f'/{sample_id}_2.10k_val_2.fq.gz'
forward_trimming_report = trim_galore_outdir + f'/{sample_id}_1.10k.fastq.gz_trimming_report.txt'
reverse_trimming_report = trim_galore_outdir + f'/{sample_id}_2.10k.fastq.gz_trimming_report.txt'

# snakemake --snakefile snakemake/1.fastq-qc.snakemake.py --cores 1 trim_galore --config sample=data/SRA/SRR6399459
rule trim_galore:
    input:
        forward_reads,
        reverse_reads
    output: 
        trimmed_forward_reads,
        trimmed_reverse_reads,
        forward_trimming_report,
        reverse_trimming_report
    shell:
        """
        trim_galore --output_dir {trim_galore_outdir} --paired {forward_reads} {reverse_reads}
        """

post_trim_galore_fastqc_outdir = config["sample"] + '/post_trim_galore_fastqc'
post_trim_forward_html_report = post_trim_galore_fastqc_outdir + "/" + sample_id + "_1.10k_val_1_fastqc.html"
post_trim_reverse_html_report = post_trim_galore_fastqc_outdir + "/" + sample_id + "_2.10k_val_2_fastqc.html"
post_trim_forward_zip = post_trim_galore_fastqc_outdir + "/" + sample_id + "_1.10k_val_1_fastqc.zip"
post_trim_reverse_zip = post_trim_galore_fastqc_outdir + "/" + sample_id + "_2.10k_val_2_fastqc.zip"

# snakemake --snakefile snakemake/1.fastq-qc.snakemake.py --cores 1 post_trim_galore_fastqc --config sample=data/SRA/SRR6399459
rule post_trim_galore_fastqc:
    input:
        trimmed_forward_reads,
        trimmed_reverse_reads
    output:
        post_trim_forward_html_report,
        post_trim_reverse_html_report,
        post_trim_forward_zip,
        post_trim_reverse_zip
    shell:
        """
        mkdir -p {post_trim_galore_fastqc_outdir}
        fastqc --outdir {post_trim_galore_fastqc_outdir} {trimmed_forward_reads} {trimmed_reverse_reads}
        """

# snakemake --snakefile snakemake/1.fastq-qc.snakemake.py --cores 1 all --config sample=data/SRA/SRR6399459
rule all:
    input:
        forward_html_report,
        reverse_html_report,
        forward_zip,
        reverse_zip,
        trimmed_forward_reads,
        trimmed_reverse_reads,
        forward_trimming_report,
        reverse_trimming_report,
        post_trim_forward_html_report,
        post_trim_reverse_html_report,
        post_trim_forward_zip,
        post_trim_reverse_zip
    shell:
        """
        echo "all done!"
        """