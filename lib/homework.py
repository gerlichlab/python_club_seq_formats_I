
def get_nucleotide_content(fastaContent):
    """Takes a fasta fileiterator created by Bio.SeqIO.parse(fastafile, "fasta")
    and returns a dictionary hast the nucleotide content of each entry
    as a fraction of total sequence lengh.
    """
    raise NotImplementedError

def remove_adapter(fasta_sequence, adapter_sequence):
    """Takes a fasta sequence (SeqRecord)
    and removes all bases after a perfect match with the adapter_sequence."""
    raise NotImplementedError