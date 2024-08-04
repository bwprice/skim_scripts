import csv
import sys
import time
import os
from Bio import Entrez, SeqIO

# Set your email and API key
Entrez.email = "example@email.com" # Add your email here
Entrez.api_key = "xxxx"  # Add your NCBI API key here

# change retmax to the number of records you want back

def fetch_refseq_protein_sequences_by_taxonomy(taxonomy, gene_name, retmax=1):
    for rank in taxonomy:
        if rank:
            search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND refseq[filter]"
            try:
                # Search for the gene in NCBI RefSeq using the constructed term
                search_handle = Entrez.esearch(db="protein", term=search_term, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    # Fetch the sequences
                    ids = search_results["IdList"]
                    print(f"Fetching sequences with IDs: {ids}")  # Debugging line
                    fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()
                    
                    # Trim the sequence header to only the accession number and add process ID as prefix
                    for record in records:
                        record.id = f"{process_id}_{record.annotations['accessions'][0]}"
                        record.description = ""
                    
                    return records, rank
                else:
                    print(f"No sequences found for {rank}. Trying next level...")
            except Exception as e:
                print(f"Error fetching data for {rank}: {e}")
            
            # Add a delay between searches if needed
            time.sleep(0)
    
    print(f"No sequences found for taxonomy {taxonomy} and gene {gene_name}.")
    return [], None

# Read the CSV file
input_file = sys.argv[1]
gene_name = sys.argv[2]  # Get the gene name from command line arguments
summary_output = []

# Create a directory for the output files
os.makedirs("protein_references", exist_ok=True)

with open(input_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        taxonomy = [row['Species'], row['Genus'], row['Family'], row['Order'], row['Class'], row['Phylum']]
        process_id = row['Process ID']
        records, matched_rank = fetch_refseq_protein_sequences_by_taxonomy(taxonomy, gene_name, retmax=1)
        
        # Save sequences to a FASTA file with the process ID and accession number as the filename
        if records:
            for record in records:
                with open(f"protein_references/{record.id}.fasta", "w") as output_handle:
                    SeqIO.write([record], output_handle, "fasta")
            
            # Add to summary output
            for record in records:
                summary_output.append({
                    'Process ID': process_id,
                    'Matched Term': matched_rank,
                    'Accession Number': record.id
                })
        
        # Add a delay between searches
        time.sleep(0)

# Write the summary output to a CSV file
output_file = input_file.replace(".csv", "_summary_output.csv")
with open(output_file, "w", newline='') as csvfile:
    fieldnames = ['Process ID', 'Matched Term', 'Accession Number']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    writer.writerows(summary_output)
