# souporcell

souporcell is a method for clustering mixed-genotype scRNAseq experiments by individual.

The inputs are just the possorted_genome_bam.bam (with index), and barcodes.tsv as output from cellranger.
souporcell is comprised of 6 steps with the first 3 using external tools and the final using the code provided here.
1. Remapping (minimap2, https://github.com/lh3/minimap2)
2. Calling candidate variants (freebayes, https://github.com/ekg/freebayes)
3. Cell allele counting (vartrix, https://github.com/10XGenomics/vartrix)
4. Clustering cells by genotype (souporcell.py)
5. Calling doublets (troublet)
6. Calling cluster genotypes and inferring amount of ambient RNA (consensus.py)

## 1. Remapping
We discuss the need for remapping in our manuscript (to be posted on biorxiv soon). We need to keep track of cell barcodes and and UMIs, so we first create a fastq with those items encoded in the readname.
> python renamer.py

