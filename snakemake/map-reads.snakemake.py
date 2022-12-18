import os

########################################################################################
# map reads with bwa
########################################################################################

# snakemake --use-conda  --snakefile snakemake/map-reads.snakemake.py --cores all bwa_index_reference_fasta --config reference_fasta=data/taxon_10239.genbank/10239.fna
rule bwa_index_reference_fasta:
    conda:
        "../environment.yml"
    input:
        config["reference_fasta"]
    output:
        config["reference_fasta"] + ".amb",
        config["reference_fasta"] + ".ann",
        config["reference_fasta"] + ".bwt",
        config["reference_fasta"] + ".pac"
    shell:
        """
        bwa index {input}
        """

# snakemake --snakefile Snakefile.py --cores all bwa_mem_align_SRR6399459_100k
rule bwa_mem_align_SRR6399459_100k:
    """
    bwa \
        mem \
        -t `nproc` \
        ref.fa \
        data/exposome/SRR6399459/SRR6399459_1.100k.fastq.gz \
        data/exposome/SRR6399459/SRR6399459_2.100k.fastq.gz \
        > data/exposome/SRR6399459/SRR6399459_2.100k.fastq.joint.fna.sam
    """

# snakemake --snakefile Snakefile.py --cores 1 samtools_filter_bwa_mem_align_SRR6399459_100k
rule samtools_filter_bwa_mem_align_SRR6399459_100k:
    """
    samtools \
        view \
        --with-header \
        --excl-flags 3844 \
        data/exposome/SRR6399459/SRR6399459_2.100k.fastq.joint.fna.sam \
        > data/exposome/SRR6399459/SRR6399459_2.100k.fastq.joint.fna.sam.filtered.sam
    """

########################################################################################
# map reads with bowtie2/hisat
########################################################################################


# snakemake --use-conda  --snakefile snakemake/map-reads.snakemake.py --cores all bowtie2_index_reference_fasta --config reference_fasta=data/taxon_10239.genbank/10239.deduped.fna
rule bowtie2_index_reference_fasta:
    conda:
        "../environment.yml"
    input:
        config["reference_fasta"]
    shell:
        """
        bowtie2-build --threads `nproc` {input} {input}
        """


# bowtie2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
# --al <path>, --al-gz <path>, --al-bz2 <path>