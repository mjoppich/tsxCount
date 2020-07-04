import os,sys
import random

with open(sys.argv[1], 'w') as fout:


	for i in range(0, int(sys.argv[2])):
		readlength = 1000

		randomreadlen = random.randint(500, 1000)
		randomreadpolyalen = random.randint(100,300)

		readseq = "".join([random.choice("ATCG") for i in range(0, randomreadlen)] + ['A' for i in range(0, randomreadpolyalen)])
		
		print("@seq{}".format(i), file=fout)
		print(readseq, file=fout)
		print("+", file=fout)
		print("".join(['&'] * len(readseq)), file=fout)

