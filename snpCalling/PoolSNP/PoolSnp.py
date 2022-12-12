import sys
from collections import defaultdict as d
import re
from optparse import OptionParser, OptionGroup
import math
import gzip
from datetime import datetime

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="""python %prog \
      --sync data.sync \
      --min-cov 10 \
      --max-cov data.cov \
      --min-count 10 \
      --min-freq 0.01 \
      --posterior-prob 0.9 \
      --miss-frac 0.1 \
      --names sample1,sample2 \
      > output.vcf
      """
parser = OptionParser(usage=usage)
helptext="""

H E L P :
_________
"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--sync", dest="m", help="A sync file")
parser.add_option("--min-cov", dest="minc", help="The minimum coverage threshold: e.g. 10")
parser.add_option("--max-cov", dest="max", help="The maximum coverage threshold percentile, e.g. 0.9")
parser.add_option("--min-count", dest="mint", help="The minimum number of counts of the alternative allele across all samples pooled")
parser.add_option("--min-freq", dest="minf", help="The minimum Frequency of the alternative allele across all samples pooled")
parser.add_option("--miss-frac", dest="mis", help="The minimum Frequency of the alternative allele across all samples pooled")
parser.add_option("--names", dest="n", help="a comma separated list of names from all samples in the sync file")
parser.add_option("--SNAPE", dest="snape", help="Set flag if input is from SNAPE",action="store_true")
parser.add_option("--posterior-prob", dest="pos", help="The threshold for the Posterior Probability calculated from SNAPE for considering a site as truly polymorphic")

parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################


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

def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string=""
    if x==".:.:.:.:.:." or x=="0:0:0:0:0:0":
        return "na"
    alleles=["A","T","C","G"]
    ah=dict(zip(alleles,map(int,x.split(":")[:4])))
    for k,v in ah.items():
        string+=v*k
    return string
covh=d(lambda:d(int))

def keywithmaxvalue(x):
    ''' This function resturns the key for the maximum value in a dictionary'''
    newhash=d(list)
    for k,v in x.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]

################################## parameters ########################################

data=options.m
minimumcount=int(options.mint)
minimumcov=int(options.minc)
minimumfreq=float(options.minf)
missfrac=float(options.mis)

############################ parse sync ###########################################

# parse sync and store alternative alleles:

# Returns a datetime object containing the local date and time
dateTimeObj = datetime.now()
print("##fileformat=VCFv4.2")
print("##fileDate="+str(dateTimeObj.day)+"/"+str(dateTimeObj.month)+"/"+str(dateTimeObj.year))
if options.snape:
    print("##Source=SNAPE")
else:
    print("##Source=PoolSnp")
print("##Parameters=<ID=MinCov,Number="+str(options.minc)+",Type=Integer,Description=\"Minimum coverage per sample\">")
print("##Parameters=<ID=MaxCov,Number="+str(options.max)+",Type=Float,Description=\"Max coverage percentile per chromosome and sample\">")
if options.snape:
    print("##Parameters=<ID=PosteriorProb,Number="+str(options.pos)+",Type=Float,Description=\"The threshold for the Posterior Probability calculated from SNAPE for considering a site as truly polymorphic\">")
else:
    print("##Parameters=<ID=MinCount,Number="+str(options.mint)+",Type=Integer,Description=\"Minimum alternative allele count across all samples pooled\">")
    print("##Parameters=<ID=MinFreq,Number="+str(options.minf)+",Type=Float,Description=\"Minimum alternative allele frequency across all samples pooled\">")
print("##Parameters=<ID=MaximumMissingFraction,Number="+str(options.mis)+",Type=Float,Description=\"Maximum fraction of samples allowed that are not fullfilling all parameters\">")
print("""##INFO=<ID=ADP,Number=1,Type=Float,Description=\"Average per-sample depth of bases\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Combined read depth across all samples\">
##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">
##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Total number of allele counts of the ALT alleles\">
##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency of the ALT alleles across all samples\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference Counts\">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Alternative Counts\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##FORMAT=<ID=FREQ,Number=1,Type=FLoat,Description=\"Variant allele frequency\">""")
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(options.n.split(",")))

for l in load_data(data):
    a=l.rstrip().split()

    CHR,POS,REF = a[:3]

    libraries=a[3:]

    ## test if missing fraction of samples smaller than threshold
    if libraries.count(".:.:.:.:.:.")/float(len(libraries)) > missfrac:
        continue

    # loop through libraries
    totalalleles=d(int)
    alleles=d(lambda:d(int))

    for j in range(len(libraries)):
        alleles[j]
        nuc = sync2string(libraries[j])

        # test if seq-string is empty
        if nuc=="na":
            continue

        # read all alleles
        for i in range(len(nuc)):

            # count alternative nucleotides
            totalalleles[nuc[i].upper()]+=1
            alleles[j][nuc[i].upper()]+=1

    ## calculate allele frequencies for individuals samples and alleles
    TAF=d(list)
    for allele,counts in totalalleles.items():
        for j in range(len(libraries)):
            if sum(alleles[j].values())!=0:
                TAF[allele].append(alleles[j][allele]/sum(alleles[j].values()))

    ## only test for min-count and min-freq if raw SYNC input and not SNAPE SYNC
    if not options.snape:

        ## test if SNPs pass minimum count / minimum frequency threshold:
        for allele,counts in totalalleles.copy().items():
            # print(allele,counts,sum(TAF[allele])/len(TAF[allele]))
            if counts<minimumcount or sum(TAF[allele])/len(TAF[allele])<minimumfreq:
                del totalalleles[allele]

    ## test if site is polymorphic
    if len(totalalleles)<2:
        #print CHR,POS,"non-poly",totalalleles
        continue

    ## create output for VCF
    ALT=[]
    ## set alternative allele order:
    for i in ["A","T","C","G"]:
        if i==REF:
            continue
        if i not in totalalleles:
            continue
        ALT.append(i)


    ## set ADP,NC,GT,AD and DP
    ADP=sum(totalalleles.values())/len(libraries)
    samplelist=[]
    co=0
    miss=0
    missN=0
    NC=0

    for j in range(len(libraries)):
        ## make empty entry if no allele counts for sample
        if libraries[j] == ".:.:.:.:.:.":
            samplelist.append("./.:.:.:.:.")
            miss+=1
            NC+=1
            continue

        alleleh = alleles[j]

        # remove alleles not counted in all samples
        for k,v in alleleh.copy().items():
            if k != REF and k not in ALT or v==0:
                del alleleh[k]
        #print(alleleh,libraries[j])
        GT,AD,RD,FREQ=[],[],0,[]
        DP=sum(alleleh.values())

        ## test if sample empty:
        if len(alleleh)==0 or sum(alleleh.values())==0:
            NC+=1
            samplelist.append("./.:.:.:.:.")
            miss+=1
            continue

        # test if population is fixed for REF allele
        if len(alleleh)==1 and REF in alleleh:
            samplelist.append("0/0:"+str(DP)+":0:"+str(DP)+":0.0")
            continue

        # test if population is fixed for ALT allele
        at=0
        if len(alleleh)==1:
            for i in range(len(ALT)):
                if ALT[i] in alleleh:
                    samplelist.append(str(i+1)+"/"+str(i+1)+":0:"+str(DP)+":"+str(DP)+":1.0")
                    at=1
                    continue
        if at==1:
            continue

        ## proceed if population not fixed
        ## set REF counts
        if REF in alleleh:
            GT.append(0)
        ## set ALT counts
        for i in range(len(ALT)):
            if ALT[i] in alleleh:
                GT.append(i+1)
                AD.append(alleleh[ALT[i]])
                RD=DP-sum(AD)
                FREQ.append(round(alleleh[ALT[i]]/float(sum(alleleh.values())),2))

        samplelist.append("/".join(map(str,GT))+":"+str(RD)+":"+",".join(map(str,AD))+":"+str(DP)+":"+",".join(map(str,FREQ)))

    ## test if missing fraction of samples smaller than threshold:
    if miss/float(len(libraries))>missfrac:
        #print CHR,POS,"missing fraction",miss/float(len(libraries))
        continue
    ADP=sum(totalalleles.values())/(len(libraries)-miss)
    DP=sum(totalalleles.values())
    AF=[]
    AC=[]
    for i in ALT:
        AF.append(str(sum(TAF[i])/len(TAF[i])))
        AC.append(str(totalalleles[i]))

    ## write output
    print(CHR+"\t"+POS+"\t.\t"+REF+"\t"+",".join(ALT)+"\t.\t.\tADP="+str(ADP)+";DP="+str(DP)+";NC="+str(NC)+";AF="+",".join(AF)+";AC="+",".join(AC)+"\tGT:RD:AD:DP:FREQ\t"+"\t".join(samplelist))
