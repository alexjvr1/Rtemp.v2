#Secondary contact across the Alps

I spoke to Ben Wielstra today (27 April 2016). He suggests that adding some of the RAD data to the phylgeography paper would be good for several reasons: 
- it would make the paper better & more interesting
- Hence I could publish in a higher impact journal
- better for my CV and for the grant applications later this year. 





Some papers to look at:

doi: 10.1111/mec.13395  : secondary contact in Salmon

Streicher et al. 2014: cytb & RAD to look at secondary contact in polytypic barking frogs. 

(i) infer phylogeographic history of CH 

(ii) identify the timing and directionality of introgression events across distinct lineages 

(iii) are the contact zones geographically in the same place for mtDNA vs RAD?

(iii) determine whether the more introgressed populations have higher genetic diversity (i.e. adaptive potential)


Aim: 

We usually assume higher genetic diversity with secondary contact, which could be associated with adaptation. 

Here we aim to find loci from secondary contact that may be important for adaptation across the Alps. 

1. determine secondary contact zone & date based on mtDNA

2. Compare with secondary contact zone 

Geographic cline analysis: 

Genomic cline analysis: 

Using the combination, 

To explore this option, I will test the ddRAD data for these popuations: 

#ddRAD dataset: 

All samples sequenced for cytb for which I have ddRAD data. = 230 individuals.

I'll create 2 datasets: 

1. CH.Phyl.230:  Full dataset from pyRAD data

2. subset.CH.Phyl:  2Mil reads per indiv (or less if less was available)
	samples with <2M reads: apla01, 

				bach04, bach08
				
				bela 06
				
				bide01, bide03
				
				bnnp01
				bnnp10
				
				buel01
				
				cava01, cava07
				
				egel01cat
				
				flue08
				
				forn01
				
				fuor04
				
				gdwe01
				
				gott05
				
				grma02, grma03, grma04, grma05
				
				grsh03
				
				gruu04
				
				hdns05
				
				jagg05
				
				mart03-05
				
				prad03, prad12
				
				rotc02
				
				rusc06
				
				sali04,06
				
				saxm01, 11
				
				seji08
				
				shwe01, 02
				
				siec08
				
				star09
				
				stir04cat
				
				tana10
				
				trim04
				
				vilt02,03
				
				wise04
				




##Dataset

This part is done on the gdc server: /gdc_home4/alexjvr/CHcomplete/outfiles_1029/Phylogeography

First I need a list of all the samples I need to keep. - I will use only the samples that I sequenced for the mtDNA. 

1. List all the individuals in the VCF file
2. make a list (PhyloNames.txt) of samples to keep

```
bcftools query -l Phylogeography/CH_6.100.vcf > Phylogeography/CHallnames.txt 
```



First I need to extract the populations from the VCF data set. But the vcf file needs to be indexed first. 

Vcftools is useful for subsetting a vcf file: 

This took a long time to run! - 32min on gdcsrv2
```
vcftools --vcf CH_6.100.vcf --keep PhyloNames.txt --recode --out CH.Phylo.vcf

Eighth Header entry should be INFO: INFO    
Keeping individuals in 'keep' list
After filtering, kept 256 out of 1029 Individuals
Outputting VCF file...
After filtering, kept 5784222 out of a possible 5784222 Sites
Run Time = 1865.00 seconds
```

I'll extract all the samples that have mtDNA sequenced. 

It turns out that I don't have RAD data for all the mtDNA sequenced samples: 256/285 samples = 89.8%

Now that I have the data, I will run the normal SNP filtering on it:*** This is before I've optimised any of the depth parameters!! 


50% genotyping rate. And MAC of 3. (across 230 indivs = 3/460 = 0.65%) - I should probably increase this! Rather use a MAF of 1%
```
vcftools --vcf CH.Phylo.vcf.recode.vcf --max-missing 0.5 --mac 3 --recode --recode-INFO-all --out s1.Phylo.RAD.vcf
```

Output: 
```
Parameters as interpreted:
	--vcf CH.Phylo.vcf.recode.vcf
	--recode-INFO-all
	--mac 3
	--max-missing 0.5
	--out s1.Phylo.RAD.vcf
	--recode

After filtering, kept 256 out of 256 Individuals
Outputting VCF file...
After filtering, kept 42827 out of a possible 5784222 Sites
Run Time = 307.00 seconds
```


And check the missingness for the individuals: 
```
vcftools --vcf s1.Phylo.RAD.vcf.recode.vcf --missing-indv

mawk '!/IN/' out.imiss | cut -f5 > totalmissing

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin( $1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
```

Output: 

![alt_txt][Fig1]
[Fig1]:https://cloud.githubusercontent.com/assets/12142475/14905048/1814df06-0d62-11e6-9b86-d63c149b44ad.png


Most of the samples have <50% missing data. 


And if I try with --max-missing 0.8 followed by the rest of the filtering 

```
vcftools --vcf CH.Phylo.vcf.recode.vcf --max-missing 0.8 --mac 3 --recode --recode-INFO-all --out s1.Phylo.RAD.vcf

Parameters as interpreted:
	--vcf CH.Phylo.vcf.recode.vcf
	--recode-INFO-all
	--mac 3
	--max-missing 0.8
	--out s1.Phylo.RAD.0.8maxmiss.vcf
	--recode

After filtering, kept 256 out of 256 Individuals
Outputting VCF file...
After filtering, kept 7710 out of a possible 5784222 Sites
Run Time = 290.00 seconds
```

![alt_txt][Fig2]
[Fig2]:https://cloud.githubusercontent.com/assets/12142475/14905392/f4dfa4c8-0d64-11e6-83a1-00c6192074b1.png

I will lose 27 individuals. Of these only wise03 is from a mixed population. (PhyloDataset.mtDNA.RAD_20160428)

I will remove these individuals. (Good news is that all the brown haplotypes are still represented. 

```
vcftools --vcf s1.Phylo.RAD.0.8maxmiss.vcf.recode.vcf --remove lowDP.indiv --recode --recode-INFO-all --out subset.imiss80

After filtering, kept 230 out of 256 Individuals
Outputting VCF file...
After filtering, kept 7710 out of a possible 7710 Sites
Run Time = 1.00 seconds
```

Rename the samples in vcf file: 

First get a list of all the samples: 
```
bcftools query -l subset.imiss80.recode.vcf

```

copy and paste this to excel. And rename accordingly (I remove the "cat" and ".fq.trim"). Nano and paste into a new file. 

Paste back: 
```
bcftools bcftools reheader subset.imiss80.recode.vcf -s CH.230.newnames.txt -o CH.230.Phylo.FINAL.vcf
```




Based on my recent checks on the pyRAD data, I should also filter all SNPs with >0.7 observed Heterozygosity. 

I will do this in R using the PLINK file. (PLINK output plink.hwe has a very strange format - multiple spaces between columns - so I couldn't figure out how to cut a specific column using linux)

In R: 
```

```
Because of a problem with ~4 of the SNPs (it printed no number under A1 allele in the plink.hwe, so couldn't get read into R), I sorted everything in excel.

There are only 15 SNPs with O.Het >0.6 (i.e. 0.19%) 

1062975:73
724541:81
1344840:45
126180:97
1391708:85
1000448:60
513575:68
472871:75
229541:28
1202593:51
371437:60
888865:72
959427:119
280931:28

Remove with plink on mac (1_Phylo/input.files/plink/)
```
nano SNPstoexclude.txt

plink --file CH.Phyl.230.7710.imiss80 --exclude SNPstoexclude.txt --recodeA --out CHPhylFINAL
```

Final dataset: 

230 individuals

7572 SNPs 

0.914192 Genotyping rate



#1. Diagnostic stats: 

1. Global pairwise Fst of mtDNA and RADdata

Use Adegenet. 

I want

1. global pairwise Fst

2. Fst between populations

Data: 


Convert the genlight object to genind
```
x.mat <- as.matrix(x) # x is a genlight object
x.mat[x.mat == 0] <- "1/1" # homozygote reference
x.mat[x.mat == 1] <- "1/2" # heterozygote
x.mat[x.mat == 2] <- "2/2" # homozygote alternate
x.gid <- df2genind(x, sep = "/", ploidy = 2)

```



#2. Population structure

###Non-geographic
	
	Model-based

	1. fastStructure
	
	2. Admixture
	
	3. Structure (not tested due to lengthy run times) 
 
	
	Model-free

	1. DAPC

	2. sNMF

	3. PCA

	Basic global population structure. 
	Assumes no LD between loci, so filter for single SNP per locus. 

###Geographic prior

	Model

	1. TESS v <3  - not tested
	
	Model-free
	
	1. sPCA

	2. TESS3

###Model comparison

Per individual: Differences in individual Q-values between clustering algorithms were assessed using the Mann-Whitney-Wilcoxon test (Wilcoxon 1945; Mann & Whitney 1947). 

Overall: Genetic structure inferred by different clustering methods were compared using a Pearson's correlation coefficient (r). 

Bonferroni correction for multiple tests was applied (Weir 1996). 


##Input file

Since most/all population structure analyses assume no linkage disequilibrium between loci, I will filter to include only single SNPs per locus. 

This is done automatically in the Structure output from PyRAD (it randomly writes a single SNP from each locus to the .str input file). But I still need to do this manually when I'm using my own vcf file. 

```
vcftools --vcf CH.230.Phylo.FINAL.vcf --thin 500 --recode --recode-INFO-all --out CH.230.Phylo.Final.thin

After filtering, kept 230 out of 230 Individuals
Outputting VCF file...
After filtering, kept 2771 out of a possible 7710 Sites
Run Time = 0.00 seconds

mv CH.230.Phylo.Final.thin.recode.vcf CH.230.2771.FINAL.thin.vcf
```

So this dramatically decreased the number of SNPs in the dataset! 

Determine the Observed Heterozygosity of this new file, and filter for Ho >0.7

1. convert to .plink

2. run hwe test

3. move to mac

4. use excel to sort columns

5. create list of SNPs to remove

6. remove with Plink

```
vcftools --vcf CH.230.2771.FINAL.thin.vcf --out CH.Phyl.230.2771 --plink
plink --file CH.Phyl.230.2771 --out CH.Phyl.230.2771.plink --recodeA


plink --file CH.Phyl.230.2771 --hardy

scp -r alexjvr@gdcsrv1.ethz.ch:/gdc_home4/alexjvr/CH.Phylogenomics/CH.Phyl.230/CH.230.2771.thinned /Users/alexjvr/2016RADAnalysis/1_Phylo/input.files/CH.230.2734.THINNED
```

This leaves only 5 SNPs with Ho > 0.7: 

```
nano SNPstoexclude.txt

	229541:28

	472871:75

	513575:68

	126180:97

	1344840:45

plink --file CH.Phyl.230.2771 --exclude SNPstoexclude.txt --recodeA --out CH.Phyl.THINNED.Final
```

Final dataset

230 individuals

2729 loci

0.914 genotyping rate



##Geographic Population Structure

###1. fastStructure

Convert input to Structure format using pgdSpider. Choose the specific fastStructure format. And change marker type to SNPs. Everything else should be left as is. All the columns are in the PLINK files. 


fastStructure can be run from the Applications folder, or specify the path in bash: 

```
python structure.py -K 4 --format=str --input=CH.230.2729 --output=CH23.2729/CH.230.2729_K4.2
```

I haven't figured out how to write a script to loop through fastStructure, but change the output file for each run. I.e. I have to manually run K 1-5 x 10 runs. The runs take just a few seconds each. 




##PCA

In R. Using the Plink input file: 




 === S4 class genlight ===
 230 genotypes,  7586 binary SNPs
 Ploidy: 2
 149815 (0.09 %) missing data
 @pop: individual membership for 230 populations
 @loc.names: labels of the SNPs
 @other: a list containing: sex  phenotype  pat  mat 


###R: PCA for CH.Phylo

```
setwd("~/2016RADAnalysis/1_Phylo/input.files/plink")

library("ade4")
library("adegenet")
library("pegas")

CH.plink <- read.PLINK("CH.Phyl.230.7710.imiss80.plink.raw")
CH.plink

indNames(CH.plink)

CH436_pop.names <- read.table("CH.Phylo.PopID.ECHN.ECHS.all.csv", header=T, quote="\"")
CH436_pop.names.factors <- as.factor(CH436_pop.names$PopID) #and convert to a factor
summary(CH436_pop.names)

pop(CH.plink) <- (CH436_pop.names.factors) #assign population names from a text file
pop(CH.plink)  ##and check that they are correct

temp <- table(unlist(other(CH.plink)))
barplot(temp,main="distribution of NoAlleles per locus",xlab="Number of Alleles",ylab="Number of sites",col=heat.colors(4))

myFreq <- glMean(CH.plink)
hist(myFreq, proba=T, col="gold", xlab="Allele Frequencies", main="Distribution of (2nd) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y,*1.8, lwd=3)

##Pop structure
pca1 <- glPca(CH.plink) ##This displays a barplot of the eigenvalues and asks the user for a number of retained principal components

scatter(pca1, posi="bottomright") #scatter plot of PC 1& 2
title("PCA of CH Rana temporaria, axes 1 and 2")

s.class(pca1$scores, pop(CH.plink), col=colors()
        [c(131,131,131,131,131,131,131,131,160,160,160,160,160,150,150,150,150,134,134,134,134,134,134,134)]) 
add.scatter.eig(pca1$eig,2,1,2) abline(h=0,v=0,col="grey")

myCol <- colorplot(pca1$scores,pca1$scores,transp=T,cex=4) abline(h=0,v=0,col="grey") 
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

#NJ tree
library(ape)
tre <- nj(dist(as.matrix(CH.plink)))
tre
plot(tre, type="unrooted", use.edge.length = TRUE,
     node.pos = NULL, show.tip.label = F, show.node.label = FALSE,
     edge.color = "black")
```


![alt_txt][Fig3]
[Fig3]:https://cloud.githubusercontent.com/assets/12142475/14938749/fc372316-0ee2-11e6-92d5-de95fb6e9264.png


The NJ tree is too messy. I'll have to think about what to do with that. 


##DAPC

Var(all) = Var(withingroups) + Var(betweengroups)

DAPC optimises Var(B) while minimising Var(A)

First convert data into Principal components, and then use a representative subset of these: 

Adegenet tutorial suggests

Use the Genlight object Plink1

```
grp <- find.clusters(Plink1, max.n.clust=100)
```
This produces a graph with the variance explained by the different PCA components: 

![alt_txt][PCA.fig]
[PCA.fig]:https://cloud.githubusercontent.com/assets/12142475/15434929/4ed4f7ca-1e6e-11e6-9e59-c750ce5a7f20.png

I'm keeping all 200 components to calculate the BIC (i.e. most likely K)

![alt_txt][BIC]
[BIC]:https://cloud.githubusercontent.com/assets/12142475/15434930/4ed63ef0-1e6e-11e6-809b-b2adae0142b1.png


Because of the distribution of PCs, I have to choose a large number of PCs to explain the data (>100). This leads to overfitting of the model. The results obtained here are meaningless. How do I deal with this problem??

I had a look at this paper:

http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1638-6

They chose 60 PCs that explains 40% of the variance in order to avoid overfitting. 


###TESS3


TESS3 uses a new method to infer ancestry: Geographically constrained least-squares estimation of ancestry coefficients.  
http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12471/epdf

K is chosen by evaluating the cross-entropy criterion for each K. This method finds the minimum number of "bits" or samples from a normal probability distribution (p) that can predict a non-normal probability distribution (q). So the smaller this number is, the better the K. 


R package has been released in devtools: 

https://github.com/cayek/TESS3/blob/master/README.md

I need to use LEA to convert my data into TESS3 format: 

For this I had to upgrade R. The following link shows how to set up R-studio to use different versions of R: 

https://support.rstudio.com/hc/en-us/articles/200486138-Using-Different-Versions-of-R

LEA is a bioconductor package. 

http://www.bioconductor.org/packages/release/bioc/html/LEA.html


In R: 

```
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

library(LEA)

setwd(/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS)
output = vcf2geno("CH.230.Phylo.FINAL.vcf")
```

The conversion removed 124 loci. Not sure why: 
```

	- number of detected individuals:	230
	- number of detected loci:		7586

For SNP info, please check ./CH.230.Phylo.FINAL.vcfsnp.

124 line(s) were removed because these are not SNPs.
Please, check ./CH.230.Phylo.FINAL.removed file, for more informations.
```

But I will work with the dataset as is. 

Now I need the .coords file, which is a file with a lat & long column for each individual (no individual names). 

list all the samples in the vcf file
```
bcftools query -l CH.230.Phylo.FINAL.vcf > CH.230.7710.names
```

```
nano CH.230.Phylo.Final.coords  ##paste all the coords into this file
```

To run TESS3: 

The executable needs to be copied to the current directory
```
cp ~/Applications/TESS3-master/build/TESS3 .

./TESS3 -x CH.230.Phylo.FINAL.geno -r CH.230.Phylo.FINAL.coords -K 1 -q K1.1.Q -g K1.1.G -f K1.1.Fst -y K1.1.sum -c 0.05
```
-I can be used to select a random subset of samples. But this full dataset ran in ~10sec, so probably not necessary. 

-y = least-squares criterion

-c = percentage of the masked genotypes. (0.05 by default). If this is set, the cross-entropy criterion is calculated. 

-i = max nr of iterations. (default = 200)


Min-entropy graph

![alt_txt][All.entropy]
[All.entropy]:https://cloud.githubusercontent.com/assets/12142475/15485793/29eaae84-20f6-11e6-860b-59a3c103a675.png


Tess3 graph for K=2 and K=3

![alt_txt][All.K2]
[All.K2]:https://cloud.githubusercontent.com/assets/12142475/15486037/a4f48d2e-20f7-11e6-83d4-1a8391ada4dc.png

![alt_txt][All.K3]
[All.K3]:https://cloud.githubusercontent.com/assets/12142475/15486038/a4fd5940-20f7-11e6-8885-df5f41892540.png


###TESS3 with a subset of the data: CHS + Brown genotypes only

First I need to select these individuals from the vcf file: 

```
vcftools --vcf results.Alldata/CH.230.Phylo.FINAL.vcf --keep samples.brownCHS --recode --recode-INFO-all --out CHS.Brown.Phylo.vcf
```

and then convert to TESS3 input using R: (Use R in command line. This is R 3.2.5)
```
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

library(LEA)

setwd(/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS)
output = vcf2geno("CHS.Brown.Phylo.vcf.recode.vcf")
```
Run TESS3 for K1-10 x5
```
cp ~/Applications/TESS3-master/build/TESS3 .

./TESS3 -x CHS.Brown.Phylo.vcf.recode.vcf.geno -r CHS.Brown.coords -K 1 -q K1.1.Q -g K1.1.G -f K1.1.Fst -y K1.1.sum -c 0.05
```

From the Min-Entropy graph, K = 2 

I interpret this as the biggest change in cross-entropy scores, as they do in this tutorial: http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf


![alt_txt][CHS.cross.entropy]
[CHS.cross.entropy]:https://cloud.githubusercontent.com/assets/12142475/15485192/c477305c-20f2-11e6-9994-e27f5b04281e.png


And the TESS3 Figure for K=2: 

![alt_txt][CHS.K2]
[CHS.K2]:https://cloud.githubusercontent.com/assets/12142475/15485345/83e9be1e-20f3-11e6-9f0b-90b4679709a2.png


Get ascii file from: 

http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp

**If the code below is used, the ascii file headers all need to be changed to caps for the script to work. 

Using the following code for the plot: 

http://membres-timc.imag.fr/Olivier.Francois/TESS_Plot.html





```
####################################
######Graph of cross-entropy scores

setwd("/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS")
library(ggplot2)

CHS.entropy <- read.csv("Cross-entropy.scores.CHS.Brown.csv")
CHS.entropy <- as.data.frame(CHS.entropy)
CHS.entropy

ggplot(CHS.entropy, aes(x=CHS.entropy$K, y=CHS.entropy$Cross.entropy)) + geom_point(shape=1) + ggtitle("Cross-entropy for CHS subset") + ylab("Cross-entropy") + xlab("K")



###Graphic display of TESS output
#########
setwd("/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS")

install.packages("fields")
install.packages("RColorBrewer")
source("MapDisplay/POPSutilities.R")

Qmatrix <- read.table("CH.230.Phylo.FINAL.3.Q")
coords <- read.table("CH.230.new.coords")
plot(coords, pch = 19, xlab = "Longitude", ylab= "Latitude")
#?map
map(add = T, boundary = T, interior = T, col = "grey80")

asc.raster=("srtm_38_02.asc")
asc.raster
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)   ##constrains the map to the raster file size
maps(matrix = Qmatrix, coords, grid, method = "max", main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude")
```














###With the old dataset (i.e. multiple SNPs per locus)

Structure is running on the GDC server (GDCsrv1 & 2) and on my computer. 

K=1-5 x 5 runs each, 50-100

Copy everything to a folder on my computer 

```
scp alexjvr@gdcsrv1.ethz.ch:/gdc_home4/alexjvr/CHcomplete/outfiles_1029/Phylogeography/StructureOnGDC/50-100/Results/* .

# my data is here on the Mac: /Users/alexjvr/2016RADAnalysis/1_Phylo/Phylo.Structure/Results

zip -r Results.zip Results
```
And upload results to Structure Harvester to check the optimal K: 

![alt_txt][Str.Harvester]
[Str.Harvester]:https://cloud.githubusercontent.com/assets/12142475/14964239/1dabe638-105c-11e6-8782-b407575131be.png

![alt_txt][Str.K23]
[Str.K23]:https://cloud.githubusercontent.com/assets/12142475/14964249/2bcb5db6-105c-11e6-8c0d-b19abf4e473b.png


###Set up hierarchical structure analysis with only CHS + Brown individuals. 

```
--vcf subset.imiss80.recode.vcf --keep CHS.Brown.names.txt --recode --out CHS.Brown.imiss80

Parameters as interpreted:
	--vcf subset.imiss80.recode.vcf
	--keep CHS.Brown.names.txt
	--out CHS.Brown.imiss80
	--recode

Keeping individuals in 'keep' list
After filtering, kept 110 out of 230 Individuals
Outputting VCF file...
After filtering, kept 7710 out of a possible 7710 Sites
Run Time = 1.00 seconds
```

Copy to mac and convert to .str format using pgdspider. 

copy back to the GDCserver and start the structure runs: 

```
structure -K 1 -o /gdc_home4/alexjvr/CHcomplete/outfiles_1029/Phylogeography/StructureOnGDC/CHS.brown.50-100/CSH.Brown.K1.run1
```

![alt_txt][CHS.deltaK]
[CHS.deltaK]:https://cloud.githubusercontent.com/assets/12142475/15059461/b5846a84-12d7-11e6-9098-983ccf04be45.png


Q by elevation

![alt_txt][CHS.elev]
[CHS.elev]:https://cloud.githubusercontent.com/assets/12142475/15059519/00855f2a-12d8-11e6-81f6-9a488943367a.png


Q by lat

![alt_txt][CHS.alt]
[CHS.alt]:https://cloud.githubusercontent.com/assets/12142475/15059528/0e509f02-12d8-11e6-9443-8f8de24a8d29.png


MAP: only the CHS/Brown samples were sequenced (i.e. for pie charts where green is also present, the bar charts represent only the RAD data of the purple/brown indivs

![alt_txt][MAP]
[MAP]:https://cloud.githubusercontent.com/assets/12142475/15059539/166453e6-12d8-11e6-9d0e-e10baa7fd4cb.png


##PCA

In R. Using the Plink input file: 




 === S4 class genlight ===
 230 genotypes,  7586 binary SNPs
 Ploidy: 2
 149815 (0.09 %) missing data
 @pop: individual membership for 230 populations
 @loc.names: labels of the SNPs
 @other: a list containing: sex  phenotype  pat  mat 


###R: PCA for CH.Phylo

```
setwd("~/2016RADAnalysis/1_Phylo/input.files/plink")

library("ade4")
library("adegenet")
library("pegas")

CH.plink <- read.PLINK("CH.Phyl.230.7710.imiss80.plink.raw")
CH.plink

indNames(CH.plink)

CH436_pop.names <- read.table("CH.Phylo.PopID.ECHN.ECHS.all.csv", header=T, quote="\"")
CH436_pop.names.factors <- as.factor(CH436_pop.names$PopID) #and convert to a factor
summary(CH436_pop.names)

pop(CH.plink) <- (CH436_pop.names.factors) #assign population names from a text file
pop(CH.plink)  ##and check that they are correct

temp <- table(unlist(other(CH.plink)))
barplot(temp,main="distribution of NoAlleles per locus",xlab="Number of Alleles",ylab="Number of sites",col=heat.colors(4))

myFreq <- glMean(CH.plink)
hist(myFreq, proba=T, col="gold", xlab="Allele Frequencies", main="Distribution of (2nd) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y,*1.8, lwd=3)

##Pop structure
pca1 <- glPca(CH.plink) ##This displays a barplot of the eigenvalues and asks the user for a number of retained principal components

scatter(pca1, posi="bottomright") #scatter plot of PC 1& 2
title("PCA of CH Rana temporaria, axes 1 and 2")

s.class(pca1$scores, pop(CH.plink), col=colors()
        [c(131,131,131,131,131,131,131,131,160,160,160,160,160,150,150,150,150,134,134,134,134,134,134,134)]) 
add.scatter.eig(pca1$eig,2,1,2) abline(h=0,v=0,col="grey")

myCol <- colorplot(pca1$scores,pca1$scores,transp=T,cex=4) abline(h=0,v=0,col="grey") 
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

#NJ tree
library(ape)
tre <- nj(dist(as.matrix(CH.plink)))
tre
plot(tre, type="unrooted", use.edge.length = TRUE,
     node.pos = NULL, show.tip.label = F, show.node.label = FALSE,
     edge.color = "black")
```


![alt_txt][Fig3]
[Fig3]:https://cloud.githubusercontent.com/assets/12142475/14938749/fc372316-0ee2-11e6-92d5-de95fb6e9264.png


The NJ tree is too messy. I'll have to think about what to do with that. 


##DAPC

Var(all) = Var(withingroups) + Var(betweengroups)

DAPC optimises Var(B) while minimising Var(A)

First convert data into Principal components, and then use a representative subset of these: 

Adegenet tutorial suggests

Use the Genlight object Plink1

```
grp <- find.clusters(Plink1, max.n.clust=100)
```
This produces a graph with the variance explained by the different PCA components: 

![alt_txt][PCA.fig]
[PCA.fig]:https://cloud.githubusercontent.com/assets/12142475/15434929/4ed4f7ca-1e6e-11e6-9e59-c750ce5a7f20.png

I'm keeping all 200 components to calculate the BIC (i.e. most likely K)

![alt_txt][BIC]
[BIC]:https://cloud.githubusercontent.com/assets/12142475/15434930/4ed63ef0-1e6e-11e6-809b-b2adae0142b1.png


Because of the distribution of PCs, I have to choose a large number of PCs to explain the data (>100). This leads to overfitting of the model. The results obtained here are meaningless. How do I deal with this problem??

I had a look at this paper:

http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1638-6

They chose 60 PCs that explains 40% of the variance in order to avoid overfitting. 


###TESS3


TESS3 uses a new method to infer ancestry: Geographically constrained least-squares estimation of ancestry coefficients.  
http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12471/epdf

K is chosen by evaluating the cross-entropy criterion for each K. This method finds the minimum number of "bits" or samples from a normal probability distribution (p) that can predict a non-normal probability distribution (q). So the smaller this number is, the better the K. 


R package has been released in devtools: 

https://github.com/cayek/TESS3/blob/master/README.md

I need to use LEA to convert my data into TESS3 format: 

For this I had to upgrade R. The following link shows how to set up R-studio to use different versions of R: 

https://support.rstudio.com/hc/en-us/articles/200486138-Using-Different-Versions-of-R

LEA is a bioconductor package. 

http://www.bioconductor.org/packages/release/bioc/html/LEA.html


In R: 

```
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

library(LEA)

setwd(/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS)
output = vcf2geno("CH.230.Phylo.FINAL.vcf")
```

The conversion removed 124 loci. Not sure why: 
```

	- number of detected individuals:	230
	- number of detected loci:		7586

For SNP info, please check ./CH.230.Phylo.FINAL.vcfsnp.

124 line(s) were removed because these are not SNPs.
Please, check ./CH.230.Phylo.FINAL.removed file, for more informations.
```

But I will work with the dataset as is. 

Now I need the .coords file, which is a file with a lat & long column for each individual (no individual names). 

list all the samples in the vcf file
```
bcftools query -l CH.230.Phylo.FINAL.vcf > CH.230.7710.names
```

```
nano CH.230.Phylo.Final.coords  ##paste all the coords into this file
```

To run TESS3: 

The executable needs to be copied to the current directory
```
cp ~/Applications/TESS3-master/build/TESS3 .

./TESS3 -x CH.230.Phylo.FINAL.geno -r CH.230.Phylo.FINAL.coords -K 1 -q K1.1.Q -g K1.1.G -f K1.1.Fst -y K1.1.sum -c 0.05
```
-I can be used to select a random subset of samples. But this full dataset ran in ~10sec, so probably not necessary. 

-y = least-squares criterion

-c = percentage of the masked genotypes. (0.05 by default). If this is set, the cross-entropy criterion is calculated. 

-i = max nr of iterations. (default = 200)


Min-entropy graph

![alt_txt][All.entropy]
[All.entropy]:https://cloud.githubusercontent.com/assets/12142475/15485793/29eaae84-20f6-11e6-860b-59a3c103a675.png


Tess3 graph for K=2 and K=3

![alt_txt][All.K2]
[All.K2]:https://cloud.githubusercontent.com/assets/12142475/15486037/a4f48d2e-20f7-11e6-83d4-1a8391ada4dc.png

![alt_txt][All.K3]
[All.K3]:https://cloud.githubusercontent.com/assets/12142475/15486038/a4fd5940-20f7-11e6-8885-df5f41892540.png


###TESS3 with a subset of the data: CHS + Brown genotypes only

First I need to select these individuals from the vcf file: 

```
vcftools --vcf results.Alldata/CH.230.Phylo.FINAL.vcf --keep samples.brownCHS --recode --recode-INFO-all --out CHS.Brown.Phylo.vcf
```

and then convert to TESS3 input using R: (Use R in command line. This is R 3.2.5)
```
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

library(LEA)

setwd(/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS)
output = vcf2geno("CHS.Brown.Phylo.vcf.recode.vcf")
```
Run TESS3 for K1-10 x5
```
cp ~/Applications/TESS3-master/build/TESS3 .

./TESS3 -x CHS.Brown.Phylo.vcf.recode.vcf.geno -r CHS.Brown.coords -K 1 -q K1.1.Q -g K1.1.G -f K1.1.Fst -y K1.1.sum -c 0.05
```

From the Min-Entropy graph, K = 2 

I interpret this as the biggest change in cross-entropy scores, as they do in this tutorial: http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf


![alt_txt][CHS.cross.entropy]
[CHS.cross.entropy]:https://cloud.githubusercontent.com/assets/12142475/15485192/c477305c-20f2-11e6-9994-e27f5b04281e.png


And the TESS3 Figure for K=2: 

![alt_txt][CHS.K2]
[CHS.K2]:https://cloud.githubusercontent.com/assets/12142475/15485345/83e9be1e-20f3-11e6-9f0b-90b4679709a2.png


Get ascii file from: 

http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp

**If the code below is used, the ascii file headers all need to be changed to caps for the script to work. 

Using the following code for the plot: 

http://membres-timc.imag.fr/Olivier.Francois/TESS_Plot.html





```
####################################
######Graph of cross-entropy scores

setwd("/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS")
library(ggplot2)

CHS.entropy <- read.csv("Cross-entropy.scores.CHS.Brown.csv")
CHS.entropy <- as.data.frame(CHS.entropy)
CHS.entropy

ggplot(CHS.entropy, aes(x=CHS.entropy$K, y=CHS.entropy$Cross.entropy)) + geom_point(shape=1) + ggtitle("Cross-entropy for CHS subset") + ylab("Cross-entropy") + xlab("K")



###Graphic display of TESS output
#########
setwd("/Users/alexjvr/2016RADAnalysis/1_Phylo/TESS")

install.packages("fields")
install.packages("RColorBrewer")
source("MapDisplay/POPSutilities.R")

Qmatrix <- read.table("CH.230.Phylo.FINAL.3.Q")
coords <- read.table("CH.230.new.coords")
plot(coords, pch = 19, xlab = "Longitude", ylab= "Latitude")
#?map
map(add = T, boundary = T, interior = T, col = "grey80")

asc.raster=("srtm_38_02.asc")
asc.raster
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster,cell_value_min=0)   ##constrains the map to the raster file size
maps(matrix = Qmatrix, coords, grid, method = "max", main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude")
```

###Hybrid Index: 

Choosing parents: 

I wanted to choose parents based on 1. mtDNA haplotype, 2. Q>0.9 based on Structure and TESS3. 

I got the following results: 



|Dataset|Parent|mtDNA+Str|+TESS3|
|:--:|:--:|:--:|:--:|
|All|ParentS|27|27|
|All|ParentN|14|13|
|CHS|ParentBrown|9|4|
|CHS|ParentCHS|3|2|


It seems like the geographic prior makes a big difference in the population assignment. 
![alt_txt][tessALL]
[tessALL]:https://cloud.githubusercontent.com/assets/12142475/15514521/4e7ddff2-219e-11e6-9eca-129932b7df26.png

![alt_txt][tessCHS]
[tessCHS]:https://cloud.githubusercontent.com/assets/12142475/15514580/959fb658-219e-11e6-8a65-032fd6c2e326.png


I have to test the following things: 

**NB reference
**Frichot et al. 2014 Fast and efficient estimation of indivdiual ancestry coefficients
**http://www.genetics.org/content/genetics/196/4/973.full.pdf


####1. Compare model-free Structure (sNMF) to TESS3
 
Using this tutorial: http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf

 
These results are much more similar to the TESS3 output. Importantly, it looks like the cross-entropy scores suggest K=3 as the most likely: 

alpha 10

![alt_txt][alpha10]
[alpha10]:https://cloud.githubusercontent.com/assets/12142475/15520537/d2572938-21bb-11e6-9398-4e69ed43b276.png

alpha 100

![alt_txt][alpha100]
[alpha100]:https://cloud.githubusercontent.com/assets/12142475/15520538/d257afd4-21bb-11e6-9d6b-fc790266298e.png




####2. In TESS3, test the effect of the alpha parameter (normalised regularisation parameter; controls the geographic regularity of the ancestry estimates). 

From Frichot et al 2014, alpha is optimised after K is chosen (from runs with alpha ~0)

So I will rerun TESS3 for K=2 for the following alpha: 

1, 10, 50, 100, 500, 1000

And compare the cross-entropy scores. Alpha minises cross-entropy scores at alpha=1


![alt_txt][alpha.opt]
[alpha.opt]:https://cloud.githubusercontent.com/assets/12142475/15520474/80296892-21bb-11e6-9abc-9819c5414ab3.png





##3. Comparison between the datasets

###1. AMOVA mtDNA vs AMOVA RAD

CHS-CHN



CHS-CHN-Brown



###sPCA





###2. Comparison between pair-wise Fst tables




##4. Contact zones

mtDNA described two contact zones, the most extensive being in the East of CH. 

###Geographic cline

A geographic cline analysis is used to describe the geographic position and extent of a contact zone. 

Here I will compare the Eastern cline for mtDNA vs RAD data

1. Data needs to be collapsed into 1D data along a straight line. If the actual contact zone is broad, or unclear, the straight-line distance can be calculated from the northern/southern/ or central point. 

Transect from GB to Zurich site. Calculate straight line distance between site and this line using "haversine" method. (www.movable-type.co.uk/scripts/latlong.html.) 

2. Model cline shape using Hzar package in R


For the Clinal analyses, I will select populations along a straight line that transects the hybrid zone from North to South. Only these populations will be used for the analyses. 

1. Use sNMF and TESS3 to calculate K & assignment probabilities

2. Compare these resutls using CLUMPP

3. Use CHN & CHS parents as indivs with mtDNA + Q>0.9 assignment to the appropriate TESS/sNMF cluster. 

4. Calculate hybrid index of individuals using Introgress in R. 

5. 

###Genomic cline

####1. Subset data from EAST populations + N + GB (parent pops)


```
bcftools query -l Phylogeography/CH_6.100.vcf > Phylogeography/CHallnames.txt 

nano EAST.names
```

```
vcftools --vcf CH.230.Phylo.FINAL.vcf --keep EAST.names --recode --recode-INFO-all --out EAST
```

and then convert to TESS3 input using R: (Use R in command line. This is R 3.2.5)
```
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

library(LEA)

setwd("/Users/alexjvr/2016RADAnalysis/1_Phylo/GenomicClines/")
output = vcf2geno("EAST.recode.vcf")
```

118 individuals

7586 individuals


####2. Use sNMF and TESS3 to calculate K & assignment probabilities (in LEA)

```
obj.east = snmf("EAST.recode.geno", K = 1:10, ploidy = 2, entropy=T, alpha=100, project="new")
plot(obj.east, col="blue4", cex=1.4, pch=19)
```
I tested alpha 0.001 (default), 1, 100, 1000. Minimal cross-entropy ~0.34, but lowest for lowers alpha. 

sNMF

![alt_txt][EAST.snmf]
[EAST.snmf]:https://cloud.githubusercontent.com/assets/12142475/15525365/92632e76-21de-11e6-856c-31a0122f57c1.png

K = 2, alpha 0.001


TESS3

```
cp ~/Applications/TESS3-master/build/TESS3 .

./TESS3 -x EAST.recode.geno -r EAST.coords -K 1 -q K1.1.Q -g K1.1.G -f K1.1.Fst -y K1.1.sum -c 0.05
```

cross-entropy scores

![alt_txt][EAST.tess]
[EAST.tess]:https://cloud.githubusercontent.com/assets/12142475/15525905/ac857f44-21e2-11e6-80c6-b4ed9edc59ad.png




####3. Use CHN & CHS parents as indivs with mtDNA + Q>0.9 assignment to the appropriate TESS/sNMF cluster. 





####4. Calculate hybrid index of individuals using Introgress in R. 





Description of the contact zones: 

1. Centrality/width

2. Geographic position



