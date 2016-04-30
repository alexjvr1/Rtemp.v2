#Secondary contact across the Alps

I spoke to Ben Wielstra today (27 April 2016). He suggests that adding some of the RAD data to the phylgeography paper would be good for several reasons: 
- it would make the paper better & more interesting
- Hence I could publish in a higher impact journal
- better for my CV and for the grant applications later this year. 


Some papers to look at:

doi: 10.1111/mec.13395  : secondary contact in Salmon

Streicher et al. 2014: cytb & RAD to look at secondary contact in polytypic barking frogs. 

(i) infer phylogenetic relationships

(ii) identify the timing and directionality of introgression events across distinct lineages 

(iii) determine whether the more introgressed populations have higher genetic diversity (i.e. adaptive potential)

(iv) 



To explore this option, I will test the ddRAD data for these popuations: 

1. Structure

2. PCA


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


##Structure

Structure is running on the server (GDCsrv1) and on my computer (to check the speed). K=2, 50-100

Started eve 29 April 2016

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



