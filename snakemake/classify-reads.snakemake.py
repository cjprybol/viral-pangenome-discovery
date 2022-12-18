import os

################################################################################################
# CLASSIFY READS WITH KRAKEN
################################################################################################

# archaea, bacteria, viral, plasmid, human1, UniVec_Core
# Standard plus protozoa & fungi
# Standard plus protozoa, fungi & plant
# capped at 8
# snakemake --snakefile Snakefile.py --cores 1 download_kraken_pluspfp8_db

# snakemake --use-conda  --snakefile snakemake/classify-reads.snakemake.py --cores 1 classify_kraken --config forward= reverse= kraken_db= out_directory=data/kraken-databases
rule classify_kraken:
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
        
########################################################################################
# classify reads with centriguge
########################################################################################
# couldn't get to work

# mkdir -p data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge
# centrifuge \
#     --threads `nproc` \
#     -x data/centrifuge-databases/nt_2018_3_3/nt \
#     -1 data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz \
#     -2 data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz \
#     --al-conc-gz data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/SRR6476469.aligned.fastq.gz \
#     -S data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/classification.txt \
#     --report-file data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/centrifuge_report.tsv

########################################################################################
# classify reads with kaiju
# https://github.com/bioinformatics-centre/kaiju
########################################################################################
# dropping because of non-inclusive database

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