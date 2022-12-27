import os

########################################################################################
# classify assembled contigs with blast
########################################################################################

# TODO
# snakemake --use-conda --snakefile snakemake/classify-contigs.snakemake.py --cores 1 download_blast_nt
rule download_blast_nt:
    conda:
        "../environment.yml"
    shell:
        """
        mkdir -p data/blastdb
        cd data/blastdb
        time update_blastdb.pl --decompress nt
        """
        
# TODO
# snakemake --use-conda --snakefile snakemake/classify-contigs.snakemake.py --cores 1 download_blast_nr
# dropped because it's too big and we have other protein databases (350-400Gb)
# rule download_blast_nr:
#     conda:
#         "../environment.yml"
#     shell:
#         """
#         mkdir -p data/blastdb
#         cd data/blastdb
#         time update_blastdb.pl --decompress nr
#         """

# update_blastdb.pl --showall pretty
# refseq_protein for protein
# nt for DNA (or env_nt for even more metagenomic projects)
# time update_blastdb.pl --decompress nt
# https://github.com/ncbi/blast_plus_docs#blast-databases
# nt

# snakemake --snakefile Snakefile.py --cores 1 blast_nt_megahit_assembled_contigs_SRR6399459_100k
rule blast_nt_megahit_assembled_contigs_SRR6399459_100k:
    # log:
    #     f"snakemake/logs/{invoked_timestamp}.log"
    conda:
        "environment.yml"
    shell:
        """
        blastn \
            -db data/blastdb/nt \
            -evalue 0.001 \
            -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles' \
            -query data/exposome/SRR6399459/megahit_100k/final.contigs.fa \
            -out data/exposome/SRR6399459/megahit_100k/final.contigs.fa.blastp.txt
        """

        
########################################################################################
# classify assembled contigs with MMSeq2
# couldn't get to work :(
########################################################################################
        
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz 87 Mb
# https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz 10 gb
# https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz 35 gb
# https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz 80 gb

# https://www.uniprot.org/help/pathway
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/README
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/READMEa
# https://github.com/soedinglab/mmseqs2/wiki#how-to-compute-the-lowest-common-ancestor-lca-of-a-given-set-of-sequences
        
    
    
# mkdir -p data/mmseq2
# mmseqs databases UniProtKB/Swiss-Prot data/mmseq2/swissprot tmp - 957Mb
# mmseqs databases UniRef50 data/mmseq2/UniRef50 tmp

# SO SLOW
# mmseqs databases UniRef90 data/mmseq2/UniRef90 tmp
# mmseqs easy-taxonomy assembled-conigs.fna data/mmseq2/swissprot alnRes tmp

# data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.faa
# data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.faa

# 2 Gb memory consumed
mmseqs easy-taxonomy data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.faa data/mmseq2/swissprot data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.faa.mmseq2.easy-taxonomy.swissprot tmp

# 2 Gb memory consumed
mmseqs easy-taxonomy data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.faa data/mmseq2/swissprot data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.faa.mmseq2.easy-taxonomy.swissprot tmp

# 23 Gb memory consumed - slow!!
# 2022-12-27T14:34:24
# 2022-12-27T16:19:08
# 2 hours
mmseqs easy-taxonomy data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.faa data/mmseq2/UniRef50 data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.faa.mmseq2.easy-taxonomy.UniRef50 tmp


mmseqs easy-taxonomy data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.faa data/mmseq2/UniRef50 data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.faa.mmseq2.easy-taxonomy.UniRef50 tmp


# don't need because using prebuilt databases??
# mmseqs createdb examples/DB.fasta targetDB
# mmseqs createtaxdb targetDB tmp
# mmseqs createindex targetDB tmp

# https://github.com/soedinglab/MMseqs2/wiki#taxonomy-report-in-kraken-or-krona-style

# mmseqs databases
# Usage: mmseqs databases <name> <o:sequenceDB> <tmpDir> [options]


# yes:
#   Name                	Type      	Taxonomy	Url
# - NR - already running through Blast
# - NR                  	Aminoacid 	     yes	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA

# maybe later:
# - NT                  	Nucleotide	       -	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
# - PDB                 	Aminoacid 	       -	https://www.rcsb.org
# - PDB70               	Profile   	       -	https://github.com/soedinglab/hh-suite
# - Pfam-A.full         	Profile   	       -	https://pfam.xfam.org
# - Pfam-A.seed         	Profile   	       -	https://pfam.xfam.org
# - CDD                   Profile                -        https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
# - eggNOG              	Profile   	       -	http://eggnog5.embl.de
# - VOGDB                 Profile                -        https://vogdb.org
# - Resfinder           	Nucleotide	       -	https://cge.cbs.dtu.dk/services/ResFinder


# annotate with prodigal

# prodigal -p single -i closed-isolate.fna -a closed-isolate.fna.prodigal.orfs.faa -d closed-isolate.fna.prodigal.orfs.fna -f gff  -o closed-isolate.fna.prodigal.orfs.gff -s closed-isolate.fna.prodigal.orfs.all

# prodigal -p meta -c -i data/SRA/SRR6399459/spades/contigs.fasta -a data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.faa -d data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.fna -f gff -o data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.gff -s data/SRA/SRR6399459/spades/contigs.fasta.prodigal.orfs.all

# pyrodigal -p meta -c -i data/SRA/SRR6399459/spades/contigs.fasta -a data/SRA/SRR6399459/spades/contigs.fasta.pyrodigal.orfs.faa -d data/SRA/SRR6399459/spades/contigs.fasta.pyrodigal.orfs.fna -f gff -o data/SRA/SRR6399459/spades/contigs.fasta.pyrodigal.orfs.gff -s data/SRA/SRR6399459/spades/contigs.fasta.pyrodigal.orfs.all

# prodigal -p meta -c -i metagenome.fna -a metagenome.fna.prodigal.orfs.faa -d metagenome.fna.prodigal.orfs.fna -f gff -o metagenome.fna.prodigal.orfs.gff -s metagenome.fna.prodigal.orfs.all

# prodigal -p meta -c -i data/SRA/SRR6399459/megahit/final.contigs.fa -a data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.faa -d data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.fna -f gff -o data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.gff -s data/SRA/SRR6399459/megahit/final.contigs.fa.prodigal.orfs.all

# pyrodigal -p meta -c -i data/SRA/SRR6399459/megahit/final.contigs.fa -a data/SRA/SRR6399459/megahit/final.contigs.fa.pyrodigal.orfs.faa -d data/SRA/SRR6399459/megahit/final.contigs.fa.pyrodigal.orfs.fna -f gff -o data/SRA/SRR6399459/megahit/final.contigs.fa.pyrodigal.orfs.gff -s data/SRA/SRR6399459/megahit/final.contigs.fa.pyrodigal.orfs.all




 -s:  Write all potential genes (with scores) to the selected file.
