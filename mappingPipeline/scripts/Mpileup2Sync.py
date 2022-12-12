import sys
from collections import defaultdict as d
import re
from optparse import OptionParser, OptionGroup
import math
import gzip
import pickle

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = """
        python %prog \
        --mpileup data.mpileup \
        --base-quality-threshold 25 \
        --minIndel 5 \
        --coding 1.8 \
        --ref reference.ref \
        --output output
        """
parser = OptionParser(usage=usage)
helptext = """

H E L P :
_________

Converts an mpileup file to a sync file and adds missing positions based on the reference genome. In addition, this script creates python object files for the chromosome-specific coverage distributions and indel positions, which occur at >= minIndel counts.

"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--mpileup", dest="m", help="A mpileup file")
parser.add_option("--ref", dest="Ref", help=" The reference genome as Python object from PickleRef.py")
parser.add_option("--base-quality-threshold", dest="b", help="The Base-quality threshold ")
parser.add_option("--coding", dest="c", help="the Illumina FASTQ quality coding")
parser.add_option("--minIndel", dest="mi", help="minimum count of Indel polymorphisms")
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

def keywithmaxvalue(d):
    ''' This function resturns the key for the maximum value in a dictionary'''
    newhash=d(list)
    for k,v in d.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]


def splitter(l, n):
    ''' This generator function returns equally sized cunks of an list'''
    #credit: Meric Lieberman, 2012
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

def extract_indel(l,TH):
    ''' This function returns an Indel from a sequence string in a pileup'''
    IndelColl=[]

    for sign in ["+","-"]:
        while sign in l:
            position = l.index(sign)
            numb =""
            i = 0
            while True:
                if l[position+1+i].isdigit():
                    numb+= l[position+1+i]
                    i+=1
                else:
                    break
            seqlength = int(numb)
            sequence = l[position:position+i+1+seqlength]
            X=l.count(sequence)
            if X>=TH:
                if sign=="-":
                    IndelColl.append([-5,seqlength+5])
                else:
                    IndelColl.append([-5,5])
            seq=l[position+1+i:position+1+i+seqlength]
            l=l.replace(sequence,"")
    # l = sequence cleaned from indels; IndelColl = dictionary with all indels; IndC = total number of indels
    return l,IndelColl


def counth2sync(x):
    ''' convert countHash to sync '''
    counts=[]
    for y in ["A","T","C","G","N","D"]:
        if y in x:
            counts.append(x[y])
        else:
            counts.append(0)
    return ":".join([str(x) for x in counts])

coverage = {}
IndelPos = {}

################################## parameters ########################################

baseqthreshold=int(options.b)
phred=float(options.c)

############################ calculate PHRED cutoff  #############################

# calculate correct PHRED score cutoff: ASCII-pc

if phred >=1.0 and phred <1.8:
    pc=64
else:
    pc=33

############################ parse FASTA ###########################################
ChrLen=d(int)
print("****** READING REF ******")
with open(options.Ref, 'rb') as reffile:
    REFID = pickle.load(reffile)
for k,v in REFID.items():
    ChrLen[k]=len(v)
print("****** READING REF DONE ******")
############################ parse MPILEUP ###########################################
print(" ")
# parse mpileup and store alternative alleles:
syncout=gzip.open(options.OUT+".sync.gz","wt")
FL=0
NUM=""
print("****** PARSING MPILEUP ******")
for l in load_data(options.m):
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
    if CHR not in coverage:
        coverage[CHR]={}
    if CHR not in IndelPos:
        IndelPos[CHR]={}
    ## test if POS = INDEX+1, i.e. the next position, otherwise fill the gaps
    if int(POS)>INDEX:
        while(INDEX<int(POS)):
            syncout.write("\t".join([CHR,str(INDEX),REFID[CHR][INDEX-1],"\t".join(["0:0:0:0:0:0"]*NUM)])+"\n")
            INDEX+=1

    # loop through libraries

    alleles=d(lambda:d(int))
    div = list(splitter(a,3))
    libraries=div[1:]

    for j in range(len(libraries)):

        alleles[j]

        if len(libraries[j])!=3:
            continue

        nuc = libraries[j][1]
        qualities = libraries[j][2]

        # test if seq-string is empty
        if nuc=="*":
            continue

        # find and remove read indices and mapping quality string
        nuc = re.sub(r'\^.',r'',nuc)
        nuc = nuc.replace('$','')
        cov = 0

        # find and remove InDels
        nuc,indel=extract_indel(nuc,int(options.mi))
        if indel!=[]:
            IndelPos[CHR][int(POS)] = indel

        # test for base quality threshold (if below: ignore nucleotide)
        #syncout.write len(nuc),len(qualities)
        nuc = "".join([nuc[x] for x in range(len(nuc)) if ord(qualities[x])-pc>=baseqthreshold])
        nuc = "".join([nuc[x] for x in range(len(nuc)) if nuc[x]!="*"])

        # store coverage distribution
        if len(nuc) not in coverage[CHR]:
            coverage[CHR][len(nuc)]=1
        else:
            coverage[CHR][len(nuc)]+=1
        if -1 not in coverage[CHR]:
            coverage[CHR][-1]=1
        else:
            coverage[CHR][-1]+=1

        # read all alleles
        for i in range(len(nuc)):

            # ignore single nucleotide deletions
            if nuc[i]=="*":
                continue
            # count nucleotides similar to reference base
            if nuc[i] =="," or nuc[i] == ".":
                alleles[j][REF]+=1
                continue
            # count alternative nucleotides
            alleles[j][nuc[i].upper()]+=1
    syncL=[]
    for k,v in sorted(alleles.items()):
        syncL.append(counth2sync(v))

    ## write output
    syncout.write(CHR+"\t"+POS+"\t"+REFID[CHR][INDEX-1]+"\t"+"\t".join(syncL)+"\n")
    INDEX+=1

## finish last chromosome
if int(POS)<ChrLen[CHR]:
    INDEX=int(POS)+1
    while(INDEX<=ChrLen[CHR]):
        print
        syncout.write("\t".join([CHR,str(INDEX),REFID[CHR][INDEX-1],"\t".join(["0:0:0:0:0:0"]*NUM)])+"\n")
        INDEX+=1

## write coverage object to file
with open(options.OUT+".cov","wb") as convout:
    pickle.dump(coverage, convout)

## write indel object to file
with open(options.OUT+".indel","wb") as indout:
    pickle.dump(IndelPos, indout)
