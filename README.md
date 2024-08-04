Hello!

Here you will find a variety of small scripts used to process museum genome skim data.

**gene_fetch.py**

This takes a spreadsheet (see: gene_fetch_example.csv) and searches refseq for a single closest reference for a specified gene.
usage = "python gene_fetch.py /PATH/TO/CSV [gene name]" 
e.g. "python gene_fetch.py /data/gene_fetch_example.csv COX1

outputs
- a folder of fasta files named based on the input name and the accession number (e.g. Sample_X_Accession_Y.fasta)
- a csv summary file with the sample name, term that matched in refseq, combined filename/sequence header

note - you'll want to check the gene names in refseq and pick the best one if there are synonyms (e.g. COI, CO1, COX1)

