from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import logging

# Set up input and output file names
fasta_file = "/Users/adams/Projects/300K/MSV000091456-SCP/fasta/UP000005640_9606.fasta"
snp_file = "/Users/adams/Projects/300K/MSV000091456-SCP/fasta/A375_mutations.csv"
output_file = "/Users/adams/Projects/300K/MSV000091456-SCP/fasta/MUT_A375_HLAI_UP000005640_9606.fasta"

# Define the number of amino acids to add before and after the SNP
num_aa = 15

# Load the FASTA file
fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Load the SNP csv
df = pd.read_csv(snp_file)
snp_df = df[df["Variant Type"] == "SNP"]
mutations_df = snp_df[snp_df["Variant Info"] == "MISSENSE"]
mutations_count = len(mutations_df["Protein Change"])
logging.info(f"There are {mutations_count} missense SNPs in the mutations df.")

ids = list(fasta_dict.keys())

for i, row in mutations_df.iterrows():
    uniprot_id = str(row["Uniprot ID"])
    for fasta_id in ids:
        if uniprot_id in fasta_id:
            mutations_df.at[i, "Fasta ID"] = fasta_id
            break

not_mapped_count = mutations_count - len(mutations_df[["Uniprot ID", "Fasta ID"]])
logging.info(f"There are {not_mapped_count} Uniprot IDs that could not be mapped to the fasta")

# Remove Proteins not found in the fasta
mutations_df = mutations_df[~mutations_df["Fasta ID"].isna()]

# Create an empty dictionary to store the SNP entries
snp_fasta_dict = {}

# Loop through the mutations and extract the mutated sequences
for index, row in mutations_df.iterrows():
    protein_id = row["Fasta ID"]
    mutation_name = row["Protein Change"]
    mutation_pos = int(row["Protein Change"][3:-1])
    protein_seq = fasta_dict[protein_id].seq
    
    # Extract the number of AA before and after the mutation
    start_pos = max(0, mutation_pos - num_aa)
    end_pos = min(len(protein_seq), mutation_pos + num_aa)
    mutated_seq = protein_seq[start_pos:end_pos]
    
  # Change the amino acid at the mutation position
    old_aa = mutation_name[2]
    new_aa = mutation_name[-1]
    mutated_seq = mutated_seq[:mutation_pos-start_pos-1] + new_aa + mutated_seq[mutation_pos-start_pos:]
    
    # Create the new entry and add it to the output file
    new_id = f"SNP_{mutation_name}_{protein_id}"
    new_seq_record = SeqIO.SeqRecord(mutated_seq, id=new_id, description="")
    snp_fasta_dict[new_id] = new_seq_record

# Add the SNP entries to the original FASTA dictionary and write to a new file
fasta_dict.update(snp_fasta_dict)

with open(output_file, "w") as output_handle:
    SeqIO.write(fasta_dict.values(), output_handle, "fasta")
