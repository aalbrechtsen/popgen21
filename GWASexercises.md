# Practical - Genome-Wide Association study (GWAS)

This exercise is about Genome-Wide Association Studies (GWAS): how to perform one and some pitfalls to look out for. 

It will be conducted from the command line using the program PLINK2. 

For a detailed description of this program click [here (original program)](http://zzz.bwh.harvard.edu/plink/) and [here (version 2 of the program)](https://www.cog-genomics.org/plink/1.9/).

The exercise will be carried out at the linux server.

## Preparation - downloading the data

First close R if you have it open (e.g. using the command "quit()"). 

Then make a folder called gwasEx for this exercise by typing
```
cd ~/
mkdir gwasEx
```
Go into this folder, download a packed version of the data that we will be analyzing today and unpack it. You can do this by typing:
```
cd gwasEx
cp /ricco/data/PhDCourse/gwas/GWASex.tar.gz .
tar -xf GWASex.tar.gz
```
List the content of the folder to see what files appeared:
```
ls
```
Your folder gwasEx should now contain a subfolder called data, containing the all files you will use in this exercise.

Some of the files are data files, but note that the folder also contain a file called "plink.plot.R"", which contain R code for plotting your results.

List the content of the "data" folder:
```
ls data/
```

## Exercise A: running your first GWAS

Briefly, the GWAS data consist of SNP genotyping data from 356 individuals some of which are have a certain disease (cases) and the rest do not (controls).

To make sure the GWAS analyses will run fast the main data file (gwa.bed) is in a binary format, which is not very reader friendly. 

However, PLINK will print summary statistics about the data (number of SNPs, number of individuals, number of cases, number of controls etc) to the screen when you run an analysis.

Also, there are two additional data files, gwa.bim and gwa.fam, which are not in binary format and which contains information about the SNPs in the data and the individuals in the data, respectively 

(you can read more about the data format in the manuals linked to above - but for now this is all you need to know).
```
wc -l data/gwa.*
```
Let's try to perform a GWAS of our data, i.e. test each SNP for association with the disease.

And let's try to do it using the simplest association test for case-control data, the allelic test, which you just performed in R in exercise 1B. 

The PLINK option "--bfile data/gwa"" will specify that the data PLINK should analyse are the files in folder called "data" with the prefix "gwa".

"—assoc" specifies that we want to use perform GWAS using the allelic test 

"—adjust"" tells PLINK to output a file that includes p-values that are adjusted for multiple testing using Bonferroni correction as well as other fancier methods.

Now perform the allelic test on all the SNPs int the dataset using PLINK by typing:
```
plink  --bfile data/gwa --logistic --adjust
```

Take a look at the text PLINK prints to your screen. Specifically, note the

 - number of SNPs

 - number of individuals

 - number of cases and controls

Next, plot the results of the GWAS by typing:
```
Rscript data/plink.plot.R plink.assoc.logistic
```
This should give you several plots.

For now just look at the Manhattan plot, which can be found in the file called plink.assoc.png. You can open it using the png-viewer eog, so by typing:
```
eog plink.assoc.logistic.png
```
A bonferroni corrected p-value threshold based on an initial p-value threshold of 0.05 is shown as a dotted line on the plot.

Explain how this threshold was reached and calculate the exact threshold using your knowledge of how many SNPs you have in your dataset (NB if you want to calculate log10 in R you can use the function log10).

Using this threshold, does any of the SNPs in your dataset seem to be associated with the disease?

Do your results seem plausible? Why/why not?


## Exercise B: checking if it went OK using QQ-plot

Now look at the QQ-plot that you already generated in exercise 2A (the file called plink.assoc.QQ.png) by typing:
```
eog plink.assoc.logistic.QQ.png
```
Here the red line is the x=y line. What does this plot suggest and why?

## Exercise 2C: doing initial QC of data part 1 (sex check)

As you can see a lot can go wrong if you do not check the quality of your data! So if you want meaningful/useful output before doing the actual testing you always have to run a lot of QC before running the association tests.

One check that is worth running is a check if the indicated genders are correct. You can check this using PLINK to calculate the inbreeding coefficient on the X chromosome under the assumption that it is an autosomal chromosome.

The reason why this is interesting is that, for technical reasons PLINK represents haploid chromosomes, such as X for males, as homozygotes. So assuming the X is an autosomal chromosome will make the males look very inbred on the X where as the woman wont (since they are diploid on the X chromosome). This means that the inbreeding coefficient estimates you get will be close to 1 for men and close to 0 for women.

This gender check can be performed in PLINK using the following command:
```
plink  --bfile data/gwa --check-sex
```
The results are in the file plink.sexcheck in which the gender reported from the person who sampled the data is in the column named PEDSEX (1 means male and 2 means female) and the inbreeding coefficient is in column named F).

Check the result by typing
```
less plink.sexcheck
```
NB you can use arrows to navigate up and down in the file and close the file viewing by typing q.
```
grep PROBLEM plink.sexcheck 
```
If you observe any problems then fix them by changing the gender in the file gwa.fam (5th colunm) (NB usually one would instead get rid of these individuals because wrong gender could indicate that the phenotypes you have do not belong to the genotyped individual. However, in this case the genders were changed on purpose for the sake of this exercises so you can safely just change them back)

## Exercise 2D: doing initial QC of data part 2 (relatedness check)

Another potential problem in association studies is spurious relatedness, where some of the individuals in the sample are closely related.

Closely related individuals can be inferred using PLINK as follows:
```
plink  --bfile data/gwa --genome
```
And you can plot the results by typing:
```
Rscript data/plink.plot.R plink.genome
```
Do that and then take a look at the result by typing
```
eog plink.genomepairwise_relatedness_1.png
```
The figure shows estimates of the relatedness for all pairs of individuals.

For each pair k1 is the proportion of the genome where the pair shares 1 of their allele identical-by-descent (IBD) and k2 is the proportion of the genome where the pair shares both their alleles IBD.
 
The expected (k1,k2) values for simple relationships are shown in the figure: MZ=monozygotic twins, PO=parent offspring, FS=full sibling, HS=half sibling, C1=first cousin (or avuncular pair), C2=cousin once removed.

Are any of the individuals in your dataset closely related?
```
awk '{if ($10>0.13) print $0}' plink.genome | grep -v nan | less -S
```
What assumption in association studies is violated when individuals are related?

And last but not least: how would you recognize if the same person is included twice (this actually happens!)

## Exercise E: doing initial QC of data part 3 (check for batch bias/non-random genotyping error)

Check if there is a batch effect/non random genotyping error by using missingness as a proxy (missingness and genotyping error are highly correlated in SNP chip data).

In other words, for each SNP test if there is a significantly different amount of missing data in the cases and controls (or equivalently if there is an association between the disease status and the missingness).

Do this with PLINK by running the following command:
```
plink --bfile data/gwa --test-missing
```
View the result in the file plink.missing, where the p-value is given in the right most coloumn or generate a plot using the Rscript data/plink.plot.R by typing this command
```
Rscript data/plink.plot.R plink.missing
```
The resulting two plots (plink.missing.png and plink.missing2.png, which you can open using the png viewer eog, so by e.g. typing eog plink.missing.png)
```
eog plink.missing.png
```
If the missingness is random the p-value should be uniformly distributed between 0 and 1.

Is this the case?

Genotyping errors are often highly correlated with missingness. How do you think this will effect your association results?

## Exercise F: doing initial QC of data part 4 (check for batch bias/non-random genotyping error again)

Principal component analysis (PCA) and a very similar methods called multidimensional scaling is also often used to reveal problems in the data.

Such analyses can be used to project all the genotype information (e.g. 500,000 marker sites) down to a low number of dimensions e.g. two.

Multidimensional scaling based on this can be performed with PLINK as follows (for all individuals except for a few which has more than 20% missingness):
```
plink --bfile data/gwa --cluster --mds-plot 2 --mind 0.2
```
Run the above command and plot the results by typing
```
Rscript data/plink.plot.R plink.mds data/gwa.fam
```
The resulting plot should now be in the file plink.mds.pdf, which you can open with the pdf viewer evince (so by typing evince plink.mds.pdf).

It shows the first two dimensions and each individual is represented by a point which is colored according to the individual's disease status.

Clustering of cases and controls is an indication of batch bias. Do you see such clustering? What else could explain this clustering?

# Exercise G: try to rerun GWAS after quality filtering SNPs

We can remove many of the error prone SNPs and individuals by removing

 - SNPs that are not in HWE

 - the rare SNPs (difficult to genotype and thus genotypes are often error prone)

 - the individuals and SNPs with lots of missing data (why?)

Let us try to do rerun an association analysis where this is done:
```
plink --bfile data/gwa --logistic --adjust --out assoc2 --hwe 0.0001 --maf 0.05 --mind 0.55 --geno 0.05
```
Plot the results using
```
Rscript data/plink.plot.R assoc2.assoc.logistic
```
How does the QQ-plot look (look at figure assoc2.assoc.QQ.png)?

What does the Manhattan plot suggest (look at figure assoc2.assoc.png)?

View the p-values adjusted for multiple testing in smarter ways than Bonferroni (e.g. FDR) using
```
less assoc2.assoc.logistic.adjusted
```
Are any of the SNPs associated when we correct for multiple testing in smarter ways?

## Exercise F: another example of a GWAS caveat

For the same individuals as above we also have another phenotype. This phenotype is strongly correlated with gender. The genotyping was done independently of this phenotype so there is no batch bias. To perform association on this phenotype type
```
plink --bfile data/gwa --logistic --pheno data/pheno3.txt --adjust --out pheno3
Rscript data/plink.plot.R pheno3.assoc.logistic
```
View the plots and results. Are any of the SNP significantly associated?
```
less pheno3.assoc.logistic.adjusted
```
Try to perform the analysis using logistic regression adjusted for sex (which is done by adding the option "—sex"):
```
plink --bfile data/gwa --logistic --out pheno3_sexAdjusted --pheno data/pheno3.txt --sex
Rscript data/plink.plot.R pheno3_sexAdjusted.assoc.logistic
```
Are there any associated SNPs according to this analysis?

Some of the probes used on the chip will hybridize with multiple loci on the genome. The associated SNPs in the previous analysis all crosshybridize with the X chromosome.

Could crosshybridization explain the difference in results from the two analyses?
