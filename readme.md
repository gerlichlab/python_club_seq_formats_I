# Python club about sequencing formats I

In this first python club about sequencing formats we will talk about formats that store nucleotide level information. These file are often text-based (due to easy of early implementation and later convention) and have a wide range of information content.

## FASTA files

The FASTA format is as easy as it gets. It consists of one or multiple read entries, each of which conforms to this format:

```
    header line ---- >READNAME ADDITIONALINFORMATION
    Sequence    ---- ACGTCCGTGAGAGAGTCGTGCACGCGACGT
```
The only restriction that it enforces is that the `READNAME` - the first set of characters after the `>` sign needs to be unique. But of course this restriction is delegated to the programs that work with FASTA files and they might or might not tell you that this happened.

FASTA files are often used in situations where quality information does not make a lot of sense, for example when describing reference genomes or primers or plasmids.