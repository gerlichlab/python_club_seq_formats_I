
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

def count_point_mutations(pysam_alignment_file):
    """Takes pysam.AlignmentFile and returns dictionary with point mutation types as keys
    and the number of the found point mutation events. E.g.
        {
            "AG": 1,
            "AT": 4,
            "AC": 5,
            "CA": 1,
            "CT": 2,
            "CG": 3,
            "TA": 1,
            "TC": 2,
            "TG": 1,
            "GA": 2,
            "GC": 3,
            "GT": 4
        }
    """
    raise NotImplementedError