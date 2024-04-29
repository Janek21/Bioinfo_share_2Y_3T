from Bio import SeqIO

record = list(SeqIO.parse("CO1A2_aln.fa", "fasta"))
count = 0
count_gaps = 0
for rec in record:
    count += len(rec)
    count_gaps += rec.count('-')

print(100*count_gaps/count)

