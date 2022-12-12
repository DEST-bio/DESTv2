import sys
from collections import defaultdict as d
import re
from optparse import OptionParser, OptionGroup
import math
import gzip
import pickle

# Based on Martin Kapun work
# Modified by Maria Bogaerts

#########################################################   HELP   #########################################################################
usage = """
        python %prog \
        --input snape_output.txt \
        --ref reference.ref \
        --output output
        """
parser = OptionParser(usage=usage)
helptext = """

H E L P :
_________

Converts SNAPE output file to a sync file and adds missing positions based on the reference genome.

"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--input", dest="IN", help="SNAPE-pooled output")
parser.add_option("--ref", dest="Ref", help=" The reference genome as Python object from PickleRef.py")
parser.add_option("--output", dest="OUT", help="output prefix")

parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin.decode('ASCII')
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y


def process_line(x):
    ''' convert SNAPE output to sync '''
    if len(x) != 11 or x[7] == "":
        return "0:0:0:0:0:0\t0:0:0:0:0:0"
    else:
        chrom = x[0]
        pos = x[1]
        ref = x[2]
        ref_alt = x[7]
        alt = ""
        if len(ref_alt) == 2:
            ref_alt_1 = ref_alt[0]
            ref_alt_2 = ref_alt[1]
            if ref_alt_1 == ref:
                alt = ref_alt_2
            elif ref_alt_2 == ref:
                alt = ref_alt_1
            else:
                print()
        elif len(ref_alt) == 1:
            if ref_alt != ref:
                alt = ref_alt
            else:
                alt = "."
        elif len(ref_alt) > 2:
            print("Found triallelic in " + str(chrom) + ":" + str(pos))
            sys.exit()
        else:
            print("ERROR")
            sys.exit()
        ref_count = x[3]
        alt_count = x[4]
        ref_qual = x[5]
        alt_qual = x[6]
        prob = x[8]
        pvalue = x[9]
        mean = (x[10]).split("\n")[0]
        a_count = 0
        t_count = 0
        c_count = 0
        g_count = 0
        if ref == "A":
            a_count = ref_count
        elif ref == "T":
            t_count = ref_count
        elif ref == "C":
            c_count = ref_count
        elif ref == "G":
            g_count = ref_count
        else:
            print("Reference is N; most frequent allele is calculated in position {}\n".format(x))
            #return "0:0:0:0:0:0\t0:0:0:0:0:0"
            # sys.exit()
        if alt == "A":
            a_count = alt_count
        elif alt == "T":
            t_count = alt_count
        elif alt == "C":
            c_count = alt_count
        elif alt == "G":
            g_count = alt_count
        #else:
        #    print("Alle position not a nucleotide base in position " + str(chrom) + " " + str(pos) + "{}\n".format(x))

        return str(a_count) + ":" + str(t_count) + ":" + str(c_count) + ":" + str(g_count) + ":0:0\t" + str(ref_count) + ":" + str(alt_count) + ":" + str(ref_qual) + ":" + str(alt_qual) + ":" + str(prob) + ":" + str(pvalue) + ":" + str(mean)


############################ parse FASTA ###########################################
ChrLen=d(int)
print("****** READING REF ******")
with open(options.Ref, 'rb') as reffile:
    REFID = pickle.load(reffile)
for k,v in REFID.items():
    ChrLen[k]=len(v)
print("****** READING REF DONE ******")
############################ parse SNAPE FILE ###########################################
print(" ")
# parse SNAPE file and store alternative alleles:
syncout=gzip.open(options.OUT+".sync.gz","wt")
FL=0
NUM=""
print("****** PARSING INPUT ******")
for l in load_data(options.IN):
    if l.rstrip()=="":
        continue
    a=l.rstrip().split("\t")
    if NUM=="":
        NUM=int(len(a)/3)-1

    ## test if first line:
    if FL==0:
        CHR=a[0]
        if CHR not in REFID:
            print(CHR+" does not match with Ref database")
            sys.exit()
        print(CHR+" started")
        #REFID[CHR]=list(REFID[CHR])
        POS=int(a[1])
        if POS>1:
            INDEX=1
            while(INDEX<POS):
                syncout.write("\t".join([CHR,str(INDEX),REFID[CHR][INDEX-1],"\t".join(["0:0:0:0:0:0"]*NUM)])+"\n")
                INDEX+=1
        else:
            INDEX=1
        FL=1

    ## test if end of chromosome:
    if CHR!=a[0]:
        ## fill up if reads not available until the last position
        if int(POS)<ChrLen[CHR]:
            INDEX=int(POS)+1
            while(INDEX<=ChrLen[CHR]):
                syncout.write("\t".join([CHR,str(INDEX),REFID[CHR][INDEX-1],"\t".join(["0:0:0:0:0:0"]*NUM)])+"\n")
                INDEX+=1
        INDEX=1
        print(a[0]+" started")

    CHR,POS,REF = a[:3]
    if CHR not in REFID:
        print(CHR+" does not match with Ref database")
        sys.exit()

    ## test if POS = INDEX+1, i.e. the next position, otherwise fill the gaps
    if int(POS)>INDEX:
        while(INDEX<int(POS)):
            syncout.write("\t".join([CHR,str(INDEX),REFID[CHR][INDEX-1],"\t".join(["0:0:0:0:0:0"]*NUM)])+"\n")
            INDEX+=1

    # loop through libraries

    syncL = process_line(a)

    ## write output
    syncout.write(CHR+"\t"+POS+"\t"+REFID[CHR][INDEX-1]+"\t"+ syncL +"\n")
    INDEX+=1

## finish last chromosome
if int(POS)<ChrLen[CHR]:
    INDEX=int(POS)+1
    while(INDEX<=ChrLen[CHR]):
        print
        syncout.write("\t".join([CHR,str(INDEX),REFID[CHR][INDEX-1],"\t".join(["0:0:0:0:0:0"]*NUM)])+"\n")
        INDEX+=1
