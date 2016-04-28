#Secondary contact across the Alps

I spoke to Ben Wielstra today (27 April 2016). He suggests that adding some of the RAD data to the phylgeography paper would be good for several reasons: 
- it would make the paper better & more interesting
- Hence I could publish in a higher impact journal
- better for my CV and for the grant applications later this year. 

To explore this option, I will test the ddRAD data for these popuations: 

1. Structure

2. PCA


First I need to extract the populations from the VCF data set: 

```
vcftools --vcf <vcf_file> --keep CEUlist.txt --out outputfile_prefix --plink
```


