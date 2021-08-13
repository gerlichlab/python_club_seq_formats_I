import itertools

def get_nucleotide_content(fastaContent):
    """Takes a fasta fileiterator created by Bio.SeqIO.parse(fastafile, "fasta")
    and returns a dictionary hast the nucleotide content of each entry
    as a fraction of total sequence lengh.
    """
    output = {}
    for seq in fastaContent:
        output[seq.name] = {}
        for nuc in ["A", "C", "G", "T"]:
            output[seq.name][nuc] = seq.seq.count(nuc)/len(seq)
    return output

def remove_adapter(fasta_sequence, adapter_sequence):
    """Takes a fasta sequence (SeqRecord)
    and removes all bases after a perfect match with the adapter_sequence."""
    occurence_index = fasta_sequence.seq.find(adapter_sequence)
    if occurence_index < 0:
        return fasta_sequence
    return fasta_sequence[:occurence_index]

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
    # create target dictionary
    output = { i + j : 0 for i,j in itertools.permutations(["A","C","G","T"], 2)}
    # accumulate reads
    for read in pysam_alignment_file:
        pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
        for pair in pairs:
            ref = pair[2].upper()
            query = read.query_sequence[pair[0]].upper()
            if ref != query:
                mut = ref + query
                output[mut] += 1
    return output