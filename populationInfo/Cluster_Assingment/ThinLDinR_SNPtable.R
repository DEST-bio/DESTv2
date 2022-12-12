#-------------------------------------------------------------------------------------
#   Made by Hansong Wang, at Hawaii Cancer Center,
# 	Recovered from http://www.gettinggeneticsdone.com/2011/03/prune-gwas-data-in-rstats.html
#	Recovered on Monday april 17, 2017
#	Recovered by JCBN
#	Select a set of SNPs > dist bp (default 100kb) apart
#   Map is a matrix with colnames "chr" and "pos".
#   Matrix MUST BE sorted by chr,pos
#   Can have more columns but the colnames "chr" and "pos" can't be changed
#   The function returns a vector of row indices corresponding
#   to the rows you should pick for your subset.

# USAGE: output_snps_to_keep <- pickSNPs(allSNPs,dist=500000)
# Then SNPs_LDthinned <- allSNPs[output_snps_to_keep, ]

#-------------------------------------------------------------------------------------
pickSNPs<-function(map, dist = 100000) { 	   
	t=as.data.frame(table(map$chr))				
	vec = map$pos
	subs = c(1,rep(NA,nrow(map)-1));  # length(subs) = nrow, but the 1st element is 1 => always select the 1st snp  
	for (k in 1:nrow(t)) { #  t: count of SNPs per chr
		if (k==1) i=1 else i=sum(t[1:(k-1),2])+1; # the 1st SNP on each ch
		subs[i] = i
		stop=sum(t[1:k,2])
		while (i<stop) {
			for (j in (i+1):stop) {
				if ((vec[j]-vec[i]) > dist) {
						#cat(i, vec[i], j, vec[j],vec[j]-vec[i], x[j],'\n'); 
						subs[j]= j; 
						i=j; 
						next; # jump out of loop
				} else if (j==stop) i=stop
			}
		}
	}
	subs[!is.na(subs)]	#  row number of selected SNPs
}