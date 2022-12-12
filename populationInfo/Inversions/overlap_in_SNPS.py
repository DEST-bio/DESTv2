import sys
from optparse import OptionParser, OptionGroup
import gzip

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --source input.vcf --target genes.sync > input.sync"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
""")
#########################################################   CODE   #########################################################################

parser.add_option("--source", dest="source", help="The source file with Chromosome and position in the first and second column")
parser.add_option("--target", dest="target", help="The target file with Chromosome and position in the first and second column")
parser.add_option("--NO", dest="NO", help="Set if NO overlap",action="store_true")


parser.add_option_group(group)
(options, args) = parser.parse_args()

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y

genehash={}

pop1=load_data(options.source)
pop2=load_data(options.target)

for l in pop1:
	if len(l.split())>1:
		genehash[l.split()[0]+"_"+l.split()[1]]=l.rstrip()

for l in pop2:
	if len(l.split())>1:
		if options.NO:
			if l.split()[0]+"_"+l.split()[1] not in genehash:
				print(l.rstrip())
		else:
			if l.split()[0]+"_"+l.split()[1] in genehash:
				print(l.rstrip())
