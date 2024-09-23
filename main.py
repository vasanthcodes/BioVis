import requests

# Amino acid table to display
amino_acid_table = {
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic Acid', 'E': 'Glutamic Acid',
    'F': 'Phenylalanine', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine', 'N': 'Asparagine',
    'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine'
}

def display_amino_acid_table():
    print("\nAmino Acid Table:")
    print(f"{'One-letter Code':<15}{'Amino Acid':<20}")
    # print("One-letter Code AminoAcid")
    print("-" * 35)
    for code, name in amino_acid_table.items():
        print(f"{code:<15}{name:<20}")
    print("-" * 35)

# Main function to display the table and take user input
if __name__ == "__main__":
    # Display the amino acid table before user input
    display_amino_acid_table()

# Dictionary of amino acid molecular weights
amino_acid_weights = {
    'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13, 'F': 165.19, 
    'G': 75.07, 'H': 155.16, 'I': 131.18, 'K': 146.19, 'L': 131.18, 
    'M': 149.21, 'N': 132.12, 'P': 115.13, 'Q': 146.15, 'R': 174.20, 
    'S': 105.09, 'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
}

# Dictionary of codons for each amino acid
codon_table = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Alanine
    'C': ['TGT', 'TGC'],  # Cysteine
    'D': ['GAT', 'GAC'],  # Aspartic Acid
    'E': ['GAA', 'GAG'],  # Glutamic Acid
    'F': ['TTT', 'TTC'],  # Phenylalanine
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],  # Glycine
    'H': ['CAT', 'CAC'],  # Histidine
    'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine
    'K': ['AAA', 'AAG'],  # Lysine
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # Leucine
    'M': ['ATG'],  # Methionine
    'N': ['AAT', 'AAC'],  # Asparagine 
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],  # Proline
    'Q': ['CAA', 'CAG'],  # Glutamine
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Threonine
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Valine
    'W': ['TGG'],  # Tryptophan
    'Y': ['TAT', 'TAC']  # Tyrosine
}

# Function to add amino acid sequence
def add_amino_acid(acid_sequence):
    return acid_sequence.upper()

# Function to calculate molecular weight of a protein
def calculate_molecular_weight(sequence):
    weight = 0.0
    for acid in sequence:
        if acid in amino_acid_weights:
            weight += amino_acid_weights[acid]
        else:
            print(f"Error: '{acid}' is not a valid amino acid.")
            return None
    return weight

# Function to display molecular weight and codons
def display_molecular_weight_and_codons(sequence):
    molecular_weight = calculate_molecular_weight(sequence)
    if molecular_weight is not None:
        print(f"\nThe molecular weight of the protein is: {molecular_weight:.2f} Da")
        print("\nAmino Acid Codon Table:")
        for acid in sequence:
            if acid in codon_table:
                print(f"Amino Acid: {acid}, Possible Codons: {', '.join(codon_table[acid])}")
            else:
                print(f"Error: '{acid}' is not a valid amino acid.")
    else:
        print("Error in molecular weight calculation.")

# Function to search UniProt by amino acid sequence
def search_uniprot_by_sequence(sequence):
    url = "https://www.uniprot.org/uniprot/"
    params = {
        'query': f'sequence:{sequence}',
        'format': 'tab',
        'columns': 'id,entry name,protein names,organism,length'
    }
    response = requests.get(url, params=params)
    if response.status_code == 200 and response.text.strip():
        print("\nProtein(s) matching the sequence found in UniProt:")
        print(response.text)
    else:
        print("No matching proteins found in UniProt for the given sequence.")

# Example usage
if __name__ == "__main__":
    sequence = input("Enter the amino acid sequence: ").strip()
    
    # Check molecular weight and codons
    print("Molecular Weight and Codons Calculation:")
    display_molecular_weight_and_codons(add_amino_acid(sequence))

    # Search protein in UniProt by sequence
    print("\nSearching for matching proteins in UniProt...")
    search_uniprot_by_sequence(sequence)
