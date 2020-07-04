from Bio import SeqIO
from collections import Counter
import sys

kmerCounts = Counter()

k=int(sys.argv[2])
seenSeqs = 0

with open(sys.argv[1], 'r') as fin:

    for record in SeqIO.parse(fin, "fastq"):
        
        seq = str(record.seq)

        #print(record.id, seq)

        if len(seq) < k:
            continue

        seenSeqs += 1

        for i in range(0, len(seq)-k):
            kseq = seq[i:i+k]

            if len(kseq) != k:
                continue

            kmerCounts[kseq] += 1


with open(sys.argv[1] + "." + str(k) + ".count", 'w') as fout:
	for kmer in kmerCounts:
		print(kmer, kmerCounts[kmer], sep="\t", file=fout)
