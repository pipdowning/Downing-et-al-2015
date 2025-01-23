# Downing-et-al-2015
Data and code for: Downing PA, Cornwallis CK, &amp; Griffin AS. 2015. Sex, long life and the evolutionary transition to cooperative breeding in birds. Proceedings of the Royal Society B, 282: 20151663


Number of supplementary items: two
1. SLC_R_Code.R
2. SLC_Tables_S1-S2.xlsx


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: SLC_R_Code.R

This R script contains all the code needed to replicate the analyses (including packages and functions).
- Data manipulation (lines 30 to 49)
- Phylogenetic heritability (lines 54 to 72)
- Question 1: Is reproduction delayed in cooperative breeders? (lines 79 to 96)
- Question 2: Do cooperative breeders live longer than non-cooperative breeders? (lines 100 to 103)
- Question 3: Does high survival make the evolution of cooperative breeding more likely? (lines 107 to 225)
- Question 4: Is long life more pronounced in promiscuous cooperative breeders? (lines 229 to 247)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: SLC_Tables_S1-S2.xlsx

This Excel document contains the following sheets:

- Table S1
	+ data on mean age at first reproduction for 39 bird species (cooperative and non-cooperative)
	+ 166 data rows
	+ column descriptions:\
		A. species = latin binomial of each species (matches the Jetz et al. nomenclature)\
		B. year = age at which reproduction first takes place\
		C. number.bred = the number of individuals that first bred at each age\
		D. cumulative.bred = the cumulative number of individuals that first bred at each age\
		E. prop.bred = cumulative the proportion of individuals breeding at each age\
		F. cooperation = was the species cooperative or non-cooperative?\
		G. study details = name of the species, reference, sample size, study duration and location

- Table S2
	+ data on annual survival and promiscuity for 238 bird species (35 cooperative and 203 non-cooperative)
	+ column descriptions:\
		A. species = latin binomial of each species (matches the Jetz et al. nomenclature)\
		B. promiscuity (%) = the percentage of broods with at least one extra-group chick\
		C. sample size = the number of broods investigate to calculate promiscuity\
		D. survival (proportion) = mean annual survival\
		E. error = statistical error for the mean annual survival estimate\
		F. error type = what measure of error was take\
		G. sample size = the number of individuals used to estimate mean annual survival\
		H. study duration (years) = the length of the annual survival study\
		I. method = the method used to estimate annual survival\
		J. location = where the annual survival study was conducted\
		K. latitude = latitude for the annual survival study location\
		L. Handbook of the Birds of the World mass (grams) = species mass taken from Handbook of the Birds of the World\
		M. breeding system = was the species cooperative or non-cooperative?\
		N. justification = reasons for breeding system classification\
		O. promiscuity reference = study from which promiscuity (%) estimates were extracted\
		P. survival reference = study from which survival (proportion) estimates were extracted\
		Q. cooperation references = study from which breeding system was determined


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
