import os

########################################################################################
# classify assembled contigs with blast
########################################################################################

# TODO
# snakemake --snakefile Snakefile.py --cores 1 download_blast_nt
rule download_blast_nt:
    conda:
        "environment.yml"
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

# snakemake --snakefile Snakefile.py --cores 1 blast_nt_megahit_assembled_contigs_SRR6399459_100k
rule blast_nt_megahit_assembled_contigs_SRR6399459_100k:
    log:
        f"snakemake/logs/{invoked_timestamp}.log"
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

# $ blastn -help
# USAGE
#   blastn [-h] [-help] [-import_search_strategy filename]
#     [-export_search_strategy filename] [-task task_name] [-db database_name]
#     [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
#     [-negative_gilist filename] [-negative_seqidlist filename]
#     [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
#     [-negative_taxidlist filename] [-entrez_query entrez_query]
#     [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
#     [-subject subject_input_file] [-subject_loc range] [-query input_file]
#     [-out output_file] [-evalue evalue] [-word_size int_value]
#     [-gapopen open_penalty] [-gapextend extend_penalty]
#     [-perc_identity float_value] [-qcov_hsp_perc float_value]
#     [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
#     [-xdrop_gap_final float_value] [-searchsp int_value]
#     [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
#     [-min_raw_gapped_score int_value] [-template_type type]
#     [-template_length int_value] [-dust DUST_options]
#     [-filtering_db filtering_database]
#     [-window_masker_taxid window_masker_taxid]
#     [-window_masker_db window_masker_db] [-soft_masking soft_masking]
#     [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
#     [-best_hit_score_edge float_value] [-subject_besthit]
#     [-window_size int_value] [-off_diagonal_range int_value]
#     [-use_index boolean] [-index_name string] [-lcase_masking]
#     [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format]
#     [-show_gis] [-num_descriptions int_value] [-num_alignments int_value]
#     [-line_length line_length] [-html] [-sorthits sort_hits]
#     [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
#     [-num_threads int_value] [-mt_mode int_value] [-remote] [-version]

#  -query <File_In>
#    Input file name
#    Default = `-'
#  -db <String>
#    BLAST database name
#     * Incompatible with:  subject, subject_loc
#  -out <File_Out, file name length < 256>
#    Output file name
#    Default = `-'
#  -evalue <Real>
#    Expectation value (E) threshold for saving hits 
#    Default = `10'

#  *** BLAST-2-Sequences options
#  -subject <File_In>
#    Subject sequence(s) to search
#     * Incompatible with:  db, gilist, seqidlist, negative_gilist,
#    negative_seqidlist, taxids, taxidlist, negative_taxids, negative_taxidlist,
#    db_soft_mask, db_hard_mask

#  *** Formatting options
#  -outfmt <String>
#    alignment view options:
#      7 = Tabular with comment lines,   
#             qseqid means Query Seq-id
#                qgi means Query GI
#               qacc means Query accesion
#            qaccver means Query accesion.version
#               qlen means Query sequence length
#             sseqid means Subject Seq-id
#          sallseqid means All subject Seq-id(s), separated by a ';'
#                sgi means Subject GI
#             sallgi means All subject GIs
#               sacc means Subject accession
#            saccver means Subject accession.version
#            sallacc means All subject accessions
#               slen means Subject sequence length
#             qstart means Start of alignment in query
#               qend means End of alignment in query
#             sstart means Start of alignment in subject
#               send means End of alignment in subject
#               qseq means Aligned part of query sequence
#               sseq means Aligned part of subject sequence
#             evalue means Expect value
#           bitscore means Bit score
#              score means Raw score
#             length means Alignment length
#             pident means Percentage of identical matches
#             nident means Number of identical matches
#           mismatch means Number of mismatches
#             staxid means Subject Taxonomy ID
#           ssciname means Subject Scientific Name
#           scomname means Subject Common Name
#         sblastname means Subject Blast Name
#          sskingdom means Subject Super Kingdom
#            staxids means unique Subject Taxonomy ID(s), separated by a ';'
#                          (in numerical order)
#          sscinames means unique Subject Scientific Name(s), separated by a ';'
#          scomnames means unique Subject Common Name(s), separated by a ';'
#         sblastnames means unique Subject Blast Name(s), separated by a ';'
#                          (in alphabetical order)
#         sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
#                          (in alphabetical order) 
#             stitle means Subject Title
#         salltitles means All Subject Title(s), separated by a '<>'
#    When not provided, the default value is:
#    'qaccver saccver pident length mismatch gapopen qstart qend sstart send
#    evalue bitscore', which is equivalent to the keyword 'std'
#    The supported format specifier for option 17 is:
#                 SQ means Include Sequence Data
#                 SR means Subject as Reference Seq
#    Default = `0'
#  -remote
#    Execute search remotely?
#     * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
#    negative_gilist, negative_seqidlist, negative_taxids, negative_taxidlist,
#    subject_loc, num_threads

########################################################################################
# classify assembled contigs with MMSeq2
# couldn't get to work :(
########################################################################################
# mmseqs databases UniProtKB/Swiss-Prot data/mmseq2/swissprot tmp

# mmseqs databases --file-allocation=none UniProtKB/Swiss-Prot data/mmseq2/swissprot tmp

# mmseqs databases UniRef90 data/mmseq2/UniRef90 tmp

# mmseqs databases UniRef90 data/mmseq2/UniRef90 tmp
# mmseqs createdb examples/DB.fasta targetDB
# mmseqs createtaxdb targetDB tmp
# mmseqs createindex targetDB tmp
# mmseqs easy-taxonomy examples/QUERY.fasta targetDB alnRes tmp


########################################################################################
# classify assembled contigs with cat
########################################################################################

# mkdir CAT
# cd CAT
# wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz

# tar -xvzf CAT_prepare_20210107.tar.gz

# CAT contigs -c {contigs fasta} -d {database folder} -t {taxonomy folder}

########################################################################################
# classify assembled contigs with SprayNPray
# https://github.com/Arkadiy-Garber/SprayNPray
########################################################################################

# wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
# git clone https://github.com/Arkadiy-Garber/SprayNPray.git
# cd SprayNPray
# bash setup.sh
# bash dlDB.sh /path/to/preferred/directory/for/db
# conda activate sprayandpray
# spray-and-pray.py -g genomeContigs.fa -out genomeContigs -ref /path/to/nr/nr.faa

########################################################################################
# classify assembled contigs with kaiju
# https://github.com/bioinformatics-centre/kaiju
########################################################################################

# mkdir data/kaijudb
# cd data/kaijudb
# kaiju-makedb -s refseq

# kaiju -z `nproc` -t nodes.dmp -f refseq/kaiju_db_refseq.fmi -i firstfile.fastq -j secondfile.fastq -o kaiju.out

# kaiju2krona -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.out.krona

# ktImportText -o kaiju.out.html kaiju.out.krona

# kaiju2table -t nodes.dmp -n names.dmp -r genus -o kaiju_summary.tsv kaiju.out [kaiju2.out, ...]

# Similarly, option -c can be used to specify the threshold by absolute read count.

# Option -u disables counting unclassified reads towards the total number of reads when calculating percentages.

# Option -p will print the full taxon path instead of just the taxon name.

# kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.names.out

# Option -u will omit unclassified reads.
# Option -p will print the full taxon path instead of just the taxon name.
# Option -r will print the path containing only to the specified ranks. For example, -r phylum,genus will append the names of phylum and genus to the end of each line.

########################################################################################
# classify assembled contigs with megan6
# https://software-ab.cs.uni-tuebingen.de/download/megan6/welcome.html
########################################################################################

########################################################################################
# classify assembled contigs with sourmash
########################################################################################

# dropping because of non-inclusive DB
# https://sourmash.readthedocs.io/en/latest/tutorials-lca.html
# curl -L -o genbank-k31.lca.json.gz https://osf.io/4f8n3/download
# sourmash sketch dna -p scaled=1000,k=31 --name-from-first some-genome.fa.gz
# sourmash lca summarize --db genbank-k31.lca.json.gz \
    # --query some-genome.fa.gz.sig
