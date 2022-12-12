import sys

### Fix FASTA header to match Repeatmasker requirements
for l in open(sys.argv[1],"r"):
	if l.startswith(">"):
		if len(l.split()[0])>50:
			print(l[:50])
		else:
			print(l.split()[0])
	else:
		print(l.rstrip())
