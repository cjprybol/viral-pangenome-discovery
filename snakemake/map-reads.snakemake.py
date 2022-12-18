import os

########################################################################################
# map reads with bwa
########################################################################################

# snakemake --snakefile Snakefile.py --cores 1 index_reference_fasta
# snakemake --use-conda  --snakefile snakemake/map-reads.snakemake.py --cores all index_reference_fasta --config reference_fasta=data/taxon_10239.genbank/10239.fna
rule bwa_index_reference_fasta:
    conda:
        "../environment.yml"
    input:
        config["reference_fasta"]
    output:
        config["reference_fasta"] + ".amb",
        config["reference_fasta"] + ".ann"
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

# Total time for call to driver() for forward index: 00:48:53
# hisat2-build -p `nproc` data/taxon_10239.genbank/joint.fna data/taxon_10239.genbank/joint.fna
# hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
# --al <path>, --al-gz <path>, --al-bz2 <path>

# --al-conc-gz /path/to/file_prefix%.aligned.fastq.gz
# --no-unal
# -p/--threads

# Total time for backward call to driver() for mirror index: 00:18:55
# bowtie2-build --threads `nproc` data/taxon_10239.genbank/joint.fna data/taxon_10239.genbank/joint.fna