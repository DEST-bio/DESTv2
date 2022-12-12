import os
import sys

snape = open(sys.argv[1], "r").readlines()
vcf_out = open((sys.argv[1]).split(".")[0] + ".vcf","w")


sample_name = (((sys.argv[1]).split("/")[-1]).split(".")[0]).split("_SNAPE")[0]


### First write header:
##fileformat=VCFv4.1
##source=snape-pooled
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference Read Counts">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternative Read Counts">
##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability">
##FORMAT=<ID=P,Number=1,Type=Float,Description="P">
##FORMAT=<ID=MEAN,Number=1,Type=Float,Description="Mean value">

print >> vcf_out, "##fileformat=VCFv4.1"
print >> vcf_out, "##source=snape-pooled"
print >> vcf_out, "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference Read Counts\">"
print >> vcf_out, "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Alternative Read Counts\">"
print >> vcf_out, "##FORMAT=<ID=RQ,Number=1,Type=Integer,Description=\"Average Quality Reference Nucleotides\">"
print >> vcf_out, "##FORMAT=<ID=AQ,Number=1,Type=Integer,Description=\"Average Quality Alternative Nucleotides\">"
print >> vcf_out, "##FORMAT=<ID=PROB,Number=1,Type=Float,Description=\"Probability\">"
print >> vcf_out, "##FORMAT=<ID=P,Number=1,Type=Float,Description=\"P\">"
print >> vcf_out, "##FORMAT=<ID=MEAN,Number=1,Type=Float,Description=\"Mean value\">"
print >> vcf_out, "#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\t" + str(sample_name)


for s in snape:
    if len(s.split("\t")) != 11:
        continue
    else:
        chrom = s.split("\t")[0]
        pos = s.split("\t")[1]
        ref = s.split("\t")[2]
        ref_alt = s.split("\t")[7]
        if len(ref_alt) == 2:
            ref_alt_1 = ref_alt[0]
            ref_alt_2 = ref_alt[1]
            if ref_alt_1 == ref:
                alt = ref_alt_2
            elif ref_alt_2 == ref:
                alt = ref_alt_1
        elif len(ref_alt) == 1:
            if ref_alt != ref:
                alt = ref_alt
            if ref_alt == ref:
                alt = "."
        elif len(ref_alt) > 2:
            print("triallelic")
        ref_count = s.split("\t")[3]
        alt_count = s.split("\t")[4]
        ref_qual = s.split("\t")[5]
        alt_qual = s.split("\t")[6]
        prob = s.split("\t")[8]
        pvalue = s.split("\t")[9]
        mean = (s.split("\t")[10]).split("\n")[0]
        print >> vcf_out, str(chrom) + "\t" + str(pos) + "\t.\t" + str(ref) + "\t" + str(alt) + "\t.\t.\t.\t" +  "RD:AD:RQ:AQ:PROB:P:MEAN" + "\t" + str(ref_count) + ":" + str(alt_count) + ":" + str(ref_qual) + ":" + str(alt_qual) + ":" + str(prob) + ":" + str(pvalue) + ":" + str(mean)

vcf_out.close()
