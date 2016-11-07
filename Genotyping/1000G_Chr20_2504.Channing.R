(under UNIX BASH Shell)
# 1000 genomes example variant data file (chromosome 20)
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# 1000 genomes phenotype data file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
# Simple but fast parser program (after compilation you'll have a program called a.out)
wget https://raw.githubusercontent.com/bwlewis/1000_genomes_examples/master/parse.c
cc -O2 parse.c
================================================================================================================================================
Under R Environment
================================================================================================================================================
library(Matrix)
p = pipe("zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  | sed /^#/d  | cut  -f '10-' | ./a.out | cut -f '1-2'")
x = read.table(p, colClasses=c("integer","integer"), fill=TRUE, row.names=NULL)

# Convert to a sparse matrix of people (rows) x variant (columns)
chr20 = sparseMatrix(i=x[,2], j=x[,1], x=1.0)

# Inspect the dimensions of this matrix
print(dim(chr20))
# [1]    2504 1812841

library(irlba)
cm = colMeans(chr20)
p = irlba(chr20, nv=3, nu=3, tol=0.1, center=cm)

library(threejs)
scatterplot3js(p$u)

# Read just the header of the chromosome file to obtain the sample identifiers
ids = readLines(pipe("zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  | sed -n /^#CHROM/p | tr '\t' '\n' | tail -n +10"))
# save(ids, file="./Ids.rda")
# Download and parse the superpopulation data for each sample, order by ids
ped = read.table(url("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"),sep="\t",header=TRUE,row.names=2)[ids,6,drop=FALSE]

# Download the subpopulation and superpopulation codes
# WARNING: These links occasionally change. Beware!
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv
# pop = read.table("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv",sep="\t",header=TRUE)
pop = read.table("./20131219.populations.tsv",sep="\t",header=TRUE)
pop = pop[1:26,]
super = pop[,3]
names(super) = pop[,2]
super = factor(super)
# The last rows of pop are summary data or non-relevant:

# Map sample sub-populations to super-populations
ped$Superpopulation = super[as.character(ped$Population)]

# Plot with colors corresponding to super populations
N = length(levels(super))
scatterplot3js(p$u, col=rainbow(N)[ped$Superpopulation], size=0.5)

## Warning: closing unused connection 6 (zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  | sed -n /^#CHROM/p | tr '    ' '
## ' | tail -n +10)


