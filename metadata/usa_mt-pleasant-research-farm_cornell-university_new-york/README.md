# Soil-Microbiome

[Comparing unamended soil with soil enriched with fresh organic matter and pyrogenic organic matter](https://www.ncbi.nlm.nih.gov/biosample?term=%22geo_loc_name=USA:%20Mt.%20Pleasant%20research%20farm,%20Cornell%20University,%20New%20York%22[attr])

Click link

1. On page, click "Send to:"
2. Under "Choose Destintation", select "File"
3. Under "Format", choose "Accessions List"
4. Click "Create File"

5. Under "Find related data", next to "Database:", select "SRA"
6. Click "Find items"

1. On page, click "Send to:"
2. Under "Choose Destintation", select "File"
3. Under "Format", choose "Summary"
4. repeat, but Under "Format", choose "Accessions List"

- we should have three files now, `biosample_result.txt`, `SraAccList.txt`, and `sra_result.csv`
- We will make a folder with a descriptive sample set ID `usa_mt-pleasant-research-farm_cornell-university_new-york`
- Then place the downloaded files into this folder
- `SraAccList.txt` has a blank line at end of file that I manually deleted with vi, but could also be handled programmatically in the future


most software I need is already available pre-installed on SCG4. Add the following to ~/.bashrc
```bash
module load htslib
module load samtools
module load sratoolkit
module load fastqc
module load jellyfish
module load spades
```

prefect downloads to ~/ncbi, and I only have 10Gb of storage there. So make an ncbi directory in $PI_HOME and then ln -s to that directly from within ~ to redirect ~/ncbi to $PI_HOME/ncbi

```bash
parallel prefetch {} :::: SraAccList.txt
```

```
mkdir 1.FASTQ
parallel fastq-dump --skip-technical --clip --read-filter pass --dumpbase --gzip --readids --split-3 --outdir 1.FASTQ {} :::: SraAccList.txt
```

fastqc [-o output dir] -d -t seqfile1 seqfile2 .. seqfileN

cutadapt?
https://github.com/FelixKrueger/TrimGalore

fastqc

merge with multiqc?

gunzip -c file.fastq.gz | jellyfish count -o file.jf -m ...
jellyfish histo -o file_jf.hist -f file.jf