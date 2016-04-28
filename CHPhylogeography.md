#Secondary contact across the Alps

I spoke to Ben Wielstra today (27 April 2016). He suggests that adding some of the RAD data to the phylgeography paper would be good for several reasons: 
- it would make the paper better & more interesting
- Hence I could publish in a higher impact journal
- better for my CV and for the grant applications later this year. 

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



First I need to extract the populations from the VCF data set: 



```
vcftools --vcf CH_6.100.vcf --keep PhyloNames.txt --out Phylo.RAD.vcf
```

Where CEUlist.txt is a list of all the that should be extracted. 

I'll extract all the samples that have mtDNA sequenced. 

It turns out that I don't have RAD data for all the mtDNA sequenced samples: 256/285 samples = 89.8%

Now that I have the data, I will run the normal SNP filtering on it: 


50% genotyping rate. And MAC of 3. 
```
vcftools --vcf Phylo.RAD.vcf --max-missing 0.5 --mac 3 --recode --recode-INFO-all --out s1.Phylo.RAD.vcf
```

Output: 
```

```


And check the missingness for the individuals: 
```
vcftools --vcf subset.g5mac3dp3.recode.vcf --missing-indv

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










