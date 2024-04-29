from Bio import SeqIO

record = list(SeqIO.parse("CO1A2_aln.fa", "fasta"))
trimmed = list(SeqIO.parse("CO1A2_trim.fa", "fasta"))

count = len(record[0])
count_trim = len(trimmed[0])

print(count-count_trim)