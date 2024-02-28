from Bio import Entrez, SeqIO

def get_cds_from_ncbi(ncbi_id):
    Entrez.email = "your_email@example.com"  # Provide your email here
    handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    cds_sequences = []
    for feature in record.features:
        if feature.type == "CDS":
            cds_sequence = feature.extract(record.seq)
            cds_sequences.append(str(cds_sequence))
    return cds_sequences

def get_codon(cds_sequence, codon_number):
    codon_start = (codon_number - 1) * 3
    codon_end = codon_start + 3
    print("CDS Sequence Length:", len(cds_sequence))  # Debugging
    print("Codon Start Index:", codon_start)  # Debugging
    print("Codon End Index:", codon_end)  # Debugging
    return cds_sequence[codon_start:codon_end]

def get_amino_acid(codon_sequence):
    codon_table = {
        "TTT": "Phe", "TTC": "Phe",
        "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
        "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
        "ATG": "Met",
        "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
        "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
        "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
        "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
        "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
        "TAT": "Tyr", "TAC": "Tyr",
        "TAA": "STOP", "TAG": "STOP", "TGA": "STOP",
        "CAT": "His", "CAC": "His",
        "CAA": "Gln", "CAG": "Gln",
        "AAT": "Asn", "AAC": "Asn",
        "AAA": "Lys", "AAG": "Lys",
        "GAT": "Asp", "GAC": "Asp",
        "GAA": "Glu", "GAG": "Glu",
        "TGT": "Cys", "TGC": "Cys",
        "TGG": "Trp",
        "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
        "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    }
    return codon_table.get(codon_sequence, "Unknown")

if __name__ == "__main__":
    ncbi_id = input("Enter the NCBI nucleotide reference sequence ID: ")
    cds_sequences = get_cds_from_ncbi(ncbi_id)
    if cds_sequences:
        cds_sequence = cds_sequences[0]  # We'll take the first CDS sequence
        print("CDS Sequence:", cds_sequence)  # Debugging
        codon_number = int(input("Enter the codon number: "))
        codon_sequence = get_codon(cds_sequence, codon_number)
        amino_acid = get_amino_acid(codon_sequence)
        print(f"The amino acid for codon {codon_sequence} is: {amino_acid}")
    else:
        print("Error: No CDS sequences found for the provided NCBI ID.")
