from Bio import SeqIO
from collections import Counter

kmerCounts = Counter()

k=14
seenSeqs = 0

with open('/mnt/d/owncloud/data/tsx/usmall_t7.fastq', 'r') as fin:

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


for kmer in kmerCounts:
    print(kmer, kmerCounts[kmer], sep="\t")
