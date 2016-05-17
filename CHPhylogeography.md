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




To explore this option, I will test the ddRAD data for these popuations: 

#1. Diagnostic stats: 

1. Global pairwise Fst of mtDNA and RADdata

```

```



#2. Population structure

1. Structure

2. PCA

3. TESS3



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


50% genotyping rate. And MAC of 3. 
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


##Structure

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

```

```



###TESS3

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

./TESS3 -x CH.230.Phylo.FINAL.geno -r CH.230.Phylo.FINAL.coords -K 3 
```
-I can be used to select a random subset of samples. But this full dataset ran in ~10sec, so probably not necessary. 


Get ascii file from: 

http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp

**If the code below is used, the ascii file headers all need to be changed to caps for the script to work. 

Using the following code for the plot: 

http://membres-timc.imag.fr/Olivier.Francois/TESS_Plot.html


```
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





##3. Comparison between the datasets

###1. AMOVA mtDNA vs AMOVA RAD

CHS-CHN



CHS-CHN-Brown






###2. Comparison between pair-wise Fst tables




##4. Contact zones

Description of the contact zones: 

1. Centrality/width

2. Geographic position



