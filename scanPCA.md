# Inference of admixture and population structure

<contents>


* PCAngsd and selection

For very resent selection we can look within closely related individuals for example with in Europeans

**data:**
 - Genotype likelihoods in Beagle format
 - ~150k random SNPs with maf > 5%
 - Four EU populations with ~100 individuals in each
 - whole genome sequencing
 - depth 2-9X (1000 genome project)

CEU | Europeans in Utah (British)
GBR | Great Britain
IBS | Iberian/Spain
TSI | Italien


First lets set the paths

<example>
# NB this must be done every time you open a new terminal
ThePath=/ricco/data/PHDCourse

## copy positions and sample information
cp $ThePath/PCangsd/data/eu1000g.sample.Info .

## load the python module
## PCAngsd
PCANGSD=`echo python3 $ThePath/prog/pcangsd/pcangsd.py`

#set path to data
EU1000=$ThePath/PCangsd/data/eu1000g.small.beagle.gz
wc eu1000g.sample.Info
N=424 #one line for header
</example>



** Explore the data
*** Take a quick look at the samples data
First try to get an overview of the dataset by copying the information file and making a summary using the following:



<example>
## view first lines of sample file
head eu1000g.sample.Info
## cut first column | sort | count
cut -f 2 -d " " eu1000g.sample.Info | sed 1d| sort | uniq -c
</example>

 - How many samples from each country?

*** Explore the input data
Now let's have a look at the GL file that you have created with ANGSD. It is a "beagle format" file called all.beagle.gz - and will be the input file to PCAangsd.
The first line in this file is a header line and after that it contains a line for each locus with GLs. By using the unix command wc we can count the number of lines in the file:

<example>
gunzip -c $EU1000 | wc -l
</example>

 - Use this to find out how many loci there are GLs for in the data set?

Next, to get an idea of what the GL file contains try from the command line to print the first 9 columns of the first 7 lines of the file:

<example>
gunzip -c $EU1000 | head -n 7 | cut -f1-9 | column -t
</example>

In general, the first three columns of a beagle file contain marker name and the two alleles, allele1 and allele2, present in the locus (in beagle A=0, C=1, G=2, T=3).
All following columns contain genotype likelihoods (three columns for each individual: first GL for homozygote for allele1,
then GL for heterozygote and then GL for homozygote for allele2). Note that the GL values sum to one per site for each individuals. This is just a normalization of the genotype likelihoods in order to avoid underflow problems in the beagle software it does not mean that they are genotype probabilities.

 - Based on this, what is the most likely genotype of Ind0 in the first locus and the locus six?


** PCAngsd

Run PCangsd with to estimate the covariance matrix while jointly estimating the individuals allele frequencies


<example>
$PCANGSD -beagle $EU1000 -o EUsmall -threads 4
</example>
The program estimates the covariance matrix that can then be used for PCA. look at the output from the program

 - The algorithm might only need to low number of PCs to estimate the allele freuqencies. How many significant PCs (see MAP test in output)?

Plot the results in R

<example>
## R
 cov <- as.matrix(read.table("EUsmall.cov"))

 e<-eigen(cov)
 ID<-read.table("eu1000g.sample.Info",head=T)
 plot(e$vectors[,1:2],col=ID$POP)

 legend("topleft",fill=1:4,levels(ID$POP))
## close R after view plot
</example>


 - Does the plot look like you expected? Which populations are close and distant to each other?


Since the European individuals in 1000G are not simple homogeneous disjoint populations it is hard to use PBS/FST or similar statistics to infer selection based on populating differences. However, PCA offers a good description of the differences between individuals which out having the define disjoint groups.

Now let try to use the PC to infer selection along the genome based on the PCA

<example>
$PCANGSD -beagle $EU1000 -o EUsmall -selection -sites_save -minMaf 0
# crate file with position and chromosome
 paste <(zcat /home/albrechtsen/embo2021/PCangsd/data/eu1000g.small.beagle.gz| cut -f 1 | sed 's/\_/\t/g' | sed 1d ) EUsmall.sites  > EUsmall.sites.info
</example>

view the SNP location info that you will need to plot the results (the third column indicate if the site is used=1 or not =0)

<example>

head EUsmall.sites.info 
</example>



plot the results of the selection scan 

<example>
library(RcppCNPy,lib="/home/albrechtsen/R/x86_64-redhat-linux-gnu-library/3.6/") # Numpy library for R

## function for QQplot
qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
legend("topleft",paste("lambda=",lambda))
}

### read in seleciton statistics (chi2 distributed)
s<-npyLoad("EUsmall.selection.npy")
## make QQ plot to QC the test statistics
qqchi(s)

# convert test statistic to p-value
pval<-1-pchisq(s,1)

## read positions (hg38)
p<-read.delim("EUsmall.sites.info",colC=c("factor","integer","integer"),head=F)

names(p)<-c("chr","pos","keep")

## make manhatten plot
plot(-log10(pval),col=p$chr[p$keep==1],xlab="Chromosomes",main="Manhatten plot")

## zoom into region
 w<-range(which(pval<1e-7)) + c(-100,100)
 keep<-w[1]:w[2]
 plot(p$pos[keep],-log10(pval[keep]),col=p$chr[keep],xlab="HG38 Position chr2")

## see the position of the most significant SNP
 p$pos[which.max(s)]
</example>

see if you can make sense of the top hit based on the genome.
 - Look in [[http://genome.ucsc.edu/cgi-bin/hgGateway][UCSC browser]]
 - Choose human GRCh38/hg38
 - search for the position of the top hit and identify the genes at that loci



* Bonus PCAngsd:  Structure from low depth sequencing data
First lets try a small data set from the low depth sequencing of the 1000 genomes  project from the following populations:

ASW     |  HapMap African ancestry individuals from SW US
CEU     | European individuals
CHB     |  Han Chinese in Beijing
JPT     | Japanese individuals
YRI     | Yoruba individuals from Nigeria
MXL | Mexican individuals from LA California

.
**data:** Genotype likelhoods in beagle format for the first 50k SNPs


To save time we have already made the input file for you for this dataset and a file with population info.
<example>

## beagle genotype likelihood file
GL1000Genomes=$ThePath/admixture/data/input.gz

## copy population information file to current folder
cp $ThePath/admixture/data/pop.info .

</example>

*** Take a quick look at the data
First try to get an overview of the dataset by copying the information file and making a summary using the following:



<example>
#copy to folder
## cut first column | sort | count
cut -f 1 -d " " pop.info | sort | uniq -c
</example>

 - Which countries are the samples from and how many samples from each?


*** Explore the input data
Now let's have a look at the GL file that you have created with ANGSD. It is a "beagle format" file called all.beagle.gz - and will be the input file to PCAangsd.
The first line in this file is a header line and after that it contains a line for each locus with GLs. By using the unix command wc we can count the number of lines in the file:

<example>
gunzip -c $GL1000Genomes | wc -l
</example>

 - Use this to find out how many loci there are GLs for in the data set?

Next, to get an idea of what the GL file contains try from the command line to print the first 9 columns of the first 7 lines of the file:

<example>
gunzip -c $GL1000Genomes | head -n 7 | cut -f1-9 | column -t
</example>

In general, the first three columns of a beagle file contain marker name and the two alleles, allele1 and allele2, present in the locus (in beagle A=0, C=1, G=2, T=3).
All following columns contain genotype likelihoods (three columns for each individual: first GL for homozygote for allele1,
then GL for heterozygote and then GL for homozygote for allele2). Note that the GL values sum to one per site for each individuals. This is just a normalization of the genotype likelihoods in order to avoid underflow problems in the beagle software it does not mean that they are genotype probabilities.

 - Based on this, what is the most likely genotype of Ind0 in the first locus and the locus six?




** PCA with admixture aware priors 

Let's try to perform PCA analysis on the same 1000 genotype genotype likelihoods 
; source /usr/local/anaconda/bin/activate PCangsd


<example>
$PCANGSD -beagle $GL1000Genomes -o input
</example>

The program estimates the covariance matrix that can then be used for PCA. look at the output from the program

 - How many significant PCs (see MAP test in output)?

Plot the results in R

<example>
#open R
pop<-read.table("pop.info")

C <- as.matrix(read.table("input.cov"))
 e <- eigen(C)
pdf("PCAngsd.pdf")
plot(e$vectors[,1:2],col=pop[,1],xlab="PC1",ylab="PC2")
legend("top",fill=1:5,levels(pop[,1]))
dev.off()
## close R
</example>

view the results with

<example>

evince PCAngsd.pdf
</example>




Compare with the estimate admixture proportions (a NGSadmix analysis)

[[../html/plots/BestK3.png][plot]]

 - In the PCA plot can you identify the Mexicans with only European ancestry?
 - What about the African American with East Asian ancestry?
 - Based on the PCA would you have the same conclusion as the admixture proportions?

Try the same analysis but without estimating individual allele frequencies. This is the same as using the first iteration of the algorithm

<example>
$PCANGSD -beagle $GL1000Genomes -o input2 -iter 0
</example>

Plot the results in R
<example>
#open R
pop<-read.table("pop.info")

C <- as.matrix(read.table("input2.cov"))
 e <- eigen(C)
pdf("PCAngsd2.pdf")
plot(e$vectors[,1:2],col=pop[,1],xlab="PC1",ylab="PC2",main="joint allele frequency")
legend("top",fill=1:5,levels(pop[,1]))
dev.off()
## close R
</example>

View plot:

<example>
evince PCAngsd2.pdf
</example>


 - Do you see any difference?
 - Would any of your conclusions change? (compared to the previous PCA plot)




Let try to use the PCA to infer admixture proportions based on the first 2 principal components. For the optimization we will use a small penalty on the admixture proportions (alpha).
<example>
 $PCANGSD -beagle $GL1000Genomes -o input -admix -admix_alpha 50
</example>

Plot the results in R
<example>
library(RcppCNPy,lib="/home/albrechtsen/R/x86_64-redhat-linux-gnu-library/3.6/") # Numpy library for R
#open R
pop<-read.table("pop.info",as.is=T)

q<-npyLoad("input.admix.3.Q.npy")

## order according to population
ord<-order(pop[,1])
barplot(t(q)[,ord],col=2:10,space=0,border=NA,xlab="Individuals",ylab="Admixture proportions")
text(tapply(1:nrow(pop),pop[ord,1],mean),-0.05,unique(pop[ord,1]),xpd=T)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)

## close R
</example>
 - how does this compare to a admixture proportion analysis  (the NGSadmix analysis above)?




** Bonus Exercise: Inbreeding in the admixed individuals

Inbreeding in admixed samples is usually not possible to estimate using standard approaches.
Let's try to estimate the inbreeding coefficient of the samples using the average allele frequency

<example>
$PCANGSD -beagle $GL1000Genomes -o IB1 -inbreedSamples -iter 1
</example>


join names and results, sort the values and look at the results
<example>
paste pop.info IB1.inbreed.samples | LC_ALL=C sort -k3g
</example>
The third column is an estimate of the inbreeding coefficient (allowing for negative)

 - Does any of the individuals look more inbreed than an offspring of a pair of first cousins  ?
 - how do you interpret negative values?
 - The results will be affected by population structure - Why?
 - see any pattern of which individuals have low (negative) and high inbreeding coefficients? - can you explain the pattern?

Now let's try to estimate the inbreeding coefficient of the samples by using the individual allele frequencies predicted by the PCA

<example>
$PCANGSD -beagle $GL1000Genomes -o IB  -inbreedSamples 
</example>


join names and results, sort the values and look at the results
<example>
paste pop.info IB.inbreed.samples | LC_ALL=C sort -k3g
</example>

 - Does any of the individual look inbreed?



* Bonus Use of NGSadmix to infer admixture proportions for numerous individuals
In this exercise we will try to use NGSadmix to analyze two different NGS datasets.

** Login to the server and set paths
First open a terminal and login to the server.
Next - before running any analyses - you need to set paths to the programs and the data you will use. Do this by pasting the following into your terminal window:

<example>
# NB this must be done every time you open a new terminal

ThePath=/home/albrechtsen/embo2021/

# Set path to ANGSD program
ANGSD=$ThePath/prog/angsd/angsd

# Set path to NGSadmix
NGSadmix=$ThePath/prog/angsd/misc/NGSadmix

# Set path to a bam file list with several bam files
BAMFOLDER=$ThePath/sfs/data/smallerbams/
</example>

** First small example
We will first try to run an NGSadmix analysis of a small dataset consisting of bam files with low depth NGS data from 30 samples: 10 from Nigeria, West Africa (YRI), 10 from Japan (JPT) and 10 with European ancestry (CEU).


CEU     | Europeans (mostly of British ancestry)
JPT     | East Asian - Japanese individuals
YRI     | West African - Nigerian Yoruba individuals
<br>




Due to computation We will use a very reduced data set:
 - Input data: bam files
 - 10 individuals from each population
 - a very reduced genome 30 x 100k random regions across the autosomes
 - Each individual is sequenced at 2-6X



**Aims**:
 - to Generate Genotype likelihood files in the beagle format
 - To infer admixture proportions for low depth sequencing data




*** Make input data using ANGSD
The input to NGSadmix is genotype likelihoods (GLs). Therefore the first step of running an NGSadmix analysis if all you have are bams files is to calculate GLs. So let's start bying doing that. First make a file that contains the paths of all the 30 bam files:

<example>
find $BAMFOLDER |  grep bam$ > all.files
</example>

To see the content of the file you made type:
<example>
cat all.files
</example>

Now calculate GLs from all the BAM files using ANGSD by running the following command in the terminal:
<example>
$ANGSD -bam all.files -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 30 -minQ 20 -minInd 25 -minMaf 0.05 -doGlf 2 -out all -P 5
</example>

NOTE that this will take a bit of time to run (a few minutes).
While waiting,  let's try to understand the above command and get some info about the data. Here is an explanation of the options used in the command:
 <example>
 -bam all.files : tells ANGSD that the bam files to calculate GL from are listed in the file all.files
 -GL 2 : tells ANGSD to use the GATK genotype likelihood model
 -doMajorMinor 1 : tells ANGSD to infer the major and minor alleles
 -doMaf 1 : tells ANGSD to calculate minor allele frequencies (needed by two of the options below: -SNP_pval and -minMaf)
 -SNP_pval 2e-6 : tells ANGSD to use a p-value threshold of 2e-6 for calling SNPs
 -minMapQ 30 : tells ANGSD what to require as minimum mapping quality (quality filter)
 -minQ 20 : tells ANGSD what to require as minimum base quality (quality filter)
 -minInd 25 : tells ANGSD to only output GLs for loci with data for at least 25 individuals for each site (quality filter)
 -minMaf 0.05 : tells ANGSD to only output GLS for loci with a minimum minor allele frequency of 0.05 (quality filter)
 -doGlf 2 : tells ANGSD to write the final genotype likelihoods into a file in beagle format
 -out all : tells ANGSD to call all output files it generate "all" followed by the appropriate suffix e.g. "all.beagle.gz"
 -P 5 : tells ANDSG use 5 threads (up to 500% CPU)
 </example>

; If ANGSD hasn't finished running yet and you are tired of waiting for it to do so, then stop it (by typing Ctrl-C) and copy the resulting files (pre-made by us) by typing:

; <example>
; cp /home/ida/teaching/popgen17/admixexercise/ANGSDoutput/all* .
; </example>

*** Explore the input data
Now let's have a look at the GL file that you have created with ANGSD. It is a "beagle format" file called all.beagle.gz - and will be the input file to NGSadmix.
The first line in this file is a header line and after that it contains a line for each locus with GLs. By using the unix command wc we can count the number of lines in the file:

<example>
gunzip -c all.beagle.gz | wc -l
</example>

 - Use this to find out how many loci there are GLs for in the data set?

Next, to get an idea of what the GL file contains try from the command line to print the first 9 columns of the first 7 lines of the file:

<example>
gunzip -c all.beagle.gz | head -n 7 | cut -f1-9 | column -t
</example>

In general, the first three columns of a beagle file contain marker name and the two alleles, allele1 and allele2, present in the locus (in beagle A=0, C=1, G=2, T=3).
All following columns contain genotype likelihoods (three columns for each individual: first GL for homozygote for allele1,
then GL for heterozygote and then GL for homozygote for allele2). Note that the GL values sum to one per site for each individuals. This is just a normalization of the genotype likelihoods in order to avoid underflow problems in the beagle software it does not mean that they are genotype probabilities.

 - Based on this, what is the most likely genotype of Ind0 in the first locus and the second locus?


*** Run an analysis of the data with NGSadmix
Now you know how the input looks. Next, let's try to perform an NGSadmix analyses of the GLs typing assuming the number of ancestral populations, K, is 3:

<example>
$NGSadmix -likes all.beagle.gz -K 3 -minMaf 0.05 -seed 1 -o all
</example>

 - While waiting for the analysis to run then make sure you understand the command. If you are in doubt seek help [[http://www.popgen.dk/software/index.php/NgsAdmix#Brief_Overview][here]]. Here you can also see what other options you have when you run an NGSadmix analyses.

*** Explore the output
The output from the analysis you just ran is three files:
 -  all.beagle.gz.log (a "log file" that summarizes the analysis run)
 -  all.beagle.gz.fopt.gz (an "fopt file", which has a line for each locus that contains an estimate of the allele frequency in each of the 3 assumed ancestral populations)
 -  all.beagle.gz.qopt.gz (a "qopt file", which has a line for each individual that contains anestimate of the individual's  ancestry proportion from each of the three assumed ancestral populations).

Let's have a look at them one at a time. First, check the log file by typing

<example>
cat all.log
</example>

 - What is the log likelihood of the estimates achieved by NGSadmix (called "best like" in the log file)?

Next, check the first line of the fopt file by typing:

<example>
zcat all.fopt.gz | head -n1
</example>

 - Based on this: what is the estimated allele frequency of the first SNP in three assumed ancestral populations?

Finally, check the first line of the qopt file and thus the estimated admixture proportions for the first individuals by typing:

<example>
head -n1 all.qopt
</example>

 - Based on this: does the individual look admixed?

You can see the ID of the first individual by getting the first line of the file you created with all your original bam files in the beginning:

<example>
head -n1 all.files
</example>

 - Based on that ID, which population does the individual come from?
 - What does this suggest about what column to look for the frequencies for that population in the qopt file?
 - Based on this and the frequency estimates for the first locus that you looked at earlier, what does NGSadmix estimate the allele frequency to be at the first locus in that population?

*** Plot the admixture proportion estimates
Finally, try to make a simple plot the estimated admixture proportions for all the individuals by opening the statistical program called R (which you do by typing "R" in the terminal and pressing enter) and then copy pasting the following code:

<example>
## open R
# Get ID and pop info for each individual
s<-strsplit(basename(scan("all.files",what="theFuck")),"\\.")
pop<-sapply(s,function(x){paste(x[5],x[1],sep="_")})

# Read inferred admixture proportions
q<-read.table("all.qopt")

# Plot them (ordered by population)
ord = order(pop)
par(mar=c(7,4,1,1))
barplot(t(q)[,ord],col=c(2,1,3),names=pop[ord],las=2,ylab="Admixture proportions",cex.names=0.75)
</example>

; If for some reason no plot pops up (technical issue) you can see the plot [[admixexercisefiles/plots/all_NGSadmix_new.pdf][here]].

Note that the order of the individuals in the plot are not the same as in the qopt file. Instead, to provide a better overview, the individuals have been ordered according to the population they are from.

 - Try to explain what the plot shows (what is on the axes, what do the colors mean and so on)
 - What does the plot suggest about whether the individuals are admixed?


NB As you could tell from the number of loci included in the analysis, the above analysis is based on data from very few loci (actually we on purpose only analyzed data from a small part of the genome to make sure the analysis ran fast). Results from an analyses of data from the entire genome can be seen [[http://popgen.dk/albrecht/phdcourse/html/plots/allWholegenome_NGSadmix.pdf][here]].

; s<-strsplit(basename(scan("/home/albrecht/public/embo2015/all.files",what="theFuck")),"\\.")
; pop<-sapply(s,function(x){paste(x[5],x[1],sep="_")})
;
; # Read inferred admixture proportions
; q<-read.table("/home/albrecht/public/embo2015/allWholegenome.beagle.gz.qopt")
;
; # Plot them (ordered by population)
; ord = order(pop)
; pdf("allWholegenome_NGSadmix.pdf")
; par(mar=c(7,3,1,1))
; barplot(t(q)[,ord],col=c(3,2,1),names=pop[ord],las=2,ylab="Admixture proportions",cex.names=0.75)
; graphics.off()

 - What does that suggest about whether the individuals are admixed?


** More realistic example
Now you know how to make input data to NGSadmix, how to run NGSadmix and what the output looks like. Let's try to look at a more realistic size dataset. More specifically let's try to run NGSadmix on data from the 1000 genomes project from the following populations:

ASW     |  HapMap African ancestry individuals from SW US
CEU     | European individuals
CHB     |  Han Chinese in Beijing
JPT     | Japanese individuals
YRI     | Yoruba individuals from Nigeria
MXL | Mexican individuals from LA California

.
**data:** Genotype likelhoods in beagle format for the first 50k SNPs


To save time we have already made the input file for you for this dataset and a file with population info.

<example>
#A file with genotype likelihoods from 100 individuals in beagle format: path:
$ThePath/admixture/data/input.gz

##A file with labels that indicate which population they are sampled from:
$ThePath/admixture/data/pop.info
</example>

*** Take a quick look at the data
First try to get an overview of the dataset by copying the information file and making a summary using the following:



<example>
#copy to folder
cp $ThePath/admixture/data/pop.info .

## cut first column | sort | count
cut -f 1 -d " " pop.info | sort | uniq -c
</example>

 - Which countries are the samples from and how many samples from each?


Now look at the genotype file input.gz. It is in the same format at the file we looked at in the previous example:

 - Use the same approach as in the previous example to find out how many loci there are genotype likelihoods from


*** Run an analysis of the data with NGSadmix
 Try to start an analysis of the data with NGSadmix with K=3 (-K 3), using 1 cpu (-P 1), using only SNPs with minor allele frequency above 0.05 (-minMaf 0.05), set the seed set to 21 (-seed 21) and set the prefix of the output files to myownoutfilesK3 (-o myownoutfilesK3).

 <example>
 $NGSadmix -likes $ThePath/admixture/data/input.gz -K 3 -P 4 -minMaf 0.05 -seed 21 -o myownoutfilesK3
 </example>

;
; Try to write the command you would run if you wanted to run an analysis of the data with NGSadmix with K=3 (-K 3), using 4 cpu (-P 4), using only SNPs with minor allele frequency above 0.05 (-minMaf 0.05), set the seed set to 21 (-seed 21) and set the prefix of the output files to myownoutfilesK3 (-o myownoutfilesK3).


Next, plot the estimated admixture proportions by running the following code in R :
<example>
## open R

## read population labels and estimated admixture proportions
pop<-read.table("pop.info",as.is=T)
q<-read.table("myownoutfilesK3.qopt")

## order according to population
ord<-order(pop[,1])
barplot(t(q)[,ord],col=2:10,space=0,border=NA,xlab="Individuals",ylab="Admixture proportions")
text(tapply(1:nrow(pop),pop[ord,1],mean),-0.05,unique(pop[ord,1]),xpd=T)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)
</example>

; Again, if for some reason no plot pops up (technical issue) you can see the plot [[URL:admixexercisefiles/plots/BestK3.png][here]].

Note that - like in the previous example - the order of the individuals in the plot are not the same as in the qopt file. Instead, to provide a better overview, the individuals have been ordered according to their population labels.

 - Try to explain what the plot shows (what is on the axes, what do the colors mean and so on)
 - What does the plot suggest in terms of population structure and admixture?

; *** Assessing convergence
;  - Try to plot the results of your own NGSadmix analysis using the same R code as above just with the path and name of the qopt file changed). Does your results differ much from the results of my best run?
;  - How much lower is the likelihood from the analysis you ran yourself? (if it is not done running then just skip this question)
; - Among my runs the difference between the likelihood from run 20 (top run) and the likelihood of the run with the 10th highest likelihood is about .02. Does it look like convergence has been reached and that the top run is (close to being) the ML solution?
; - Among the runs I ran I got this result: [[admixexercisefiles/NGSadmixoutput/myoutfilesxxK3.log][log file]] and [[admixexercisefiles/NGSadmixoutput/myoutfilesxxK3.qopt][qopt file]]. Try to look at the likelihood and compare it to the likelihood from the top run. What do you think happened?
; - Try to plot the estimated admixture proportions for this run (the qopt file is located the same place as the qopt file for the best run and is myoutfilesxxK3.qopt so you can plot it by using the same R code as before but just with file name changed from myoutfiles20K3.qopt to myoutfilesxxK3.qopt).
; Again, if for some reason no plot pops up (technical issue) you can see the plot [[URL:admixexercisefiles/plots/badK3.png][here]].
; - What conclusion would you have made for the Europeans if you had only done this single run (and knew nothing about the demographic history of the populations you were analyzing)?

; *** Look at your own results (if you have time)
; - What is the likelihood for you own run? (look in the log file)?
; - Is it close to that of the best of my runs?
; - Try to plot the estimated admixture proportions
; - Is the result similar to the the best of my solutions?

*** Other K values (if you have time)
 - Try to run NGSadmix with K=4 instead.
 - Plot the output (if you have trouble plotting it here is the plot I got: [[URL:http://popgen.dk/albrecht/phdcourse/admixture/data/bestK4.png][K4 plot]].
; - What do they suggest in terms of population structure and admixture?
 - Based on all the results what can you say about the Mexican samples (MXL)?

; *** If you do not have NGS data (tip)
; If you do not have NGS data but instead called genotypes e.g. from a SNP chip another easy to use tool is Admixture.
; Another option is STRUCTURE which you can also use for microsatellite data.
; You can read about these programs here: [[http://www.genetics.ucla.edu/software/admixture/index.html][Admixture]] and [[http://pritchardlab.stanford.edu/structure.html][STRUCTURE]].
