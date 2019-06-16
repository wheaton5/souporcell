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

Easy Installation (Linux) (recommended) 
2 Options: 
1. download singularity vm image (1gb) (recommended)
2. build vm image (30min, requires root)

1. download from google drive. Google drive makes it annoyingly difficult to download via the terminal.
```

```
check if your server has singularity, most science HPC clusters should have this. Email your admins if you don't have it and you don't have root to install it.
```
which singularity
```


Mac Installation (work in progress, not currently supported)
You will need homebrew to install vagrant and vagrant to start a vm with singularity
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install vagrant
vagrant init sylabs/singularity-3.2-ubuntu-bionic64
vagrant up
```


Or you can install everything independently (not recommended, but shouldn't be too bad)
```
git clone https://github.com/wheaton5/souporcell.git
```

## 1. Remapping
We discuss the need for remapping in our manuscript (to be posted on biorxiv soon). We need to keep track of cell barcodes and and UMIs, so we first create a fastq with those items encoded in the readname.
Requires python 3.0, modules pysam, argparse (pip install/conda install depending on environment)
Easiest to first add the souporcell directory to your PATH variable with 
```
export PATH=/path/to/souporcell:$PATH
```
Then run the renamer.py script to put some of the quality information in the read name. For human data this step will typically take several hours and the output fq file will be somewhat larger than the input bam
```
python renamer.py --bam possorted_genome_bam.bam --barcodes barcodes.tsv --out fq.fq
```
Then we must remap these reads using minimap2 (similar results have been seen with hisat2)
Requires minimap2 https://github.com/lh3/minimap2
and add /path/to/minimap2 to your PATH. For human data the remapping will typically require more than 12 Gb memory and may take a few hours to run.
```
minimap2 -ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no <reference_fasta_file> <fastq_file> > minimap.sam
```
(note the -t 8 as the number of threads, change this as needed)
Now we must retag the reads with their cell barcodes and UMIs
```
python retag.py --sam minimap.sam --out minitagged.bam
```
Then we must sort and index our bam.
Requires samtools http://www.htslib.org/
```
samtools sort minitagged.bam > minitagged_sorted.bam
samtools index minitagged_sorted.bam
```

## 2. Calling candidate variants
You may wish to break this into multiple jobs such as 1 job per chromosome and merge after but the basic command is the following.
Requires freebayes https://github.com/ekg/freebayes and add /path/to/freebayes/bin to your PATH
```
freebayes -f <reference_fasta> -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --max-coverage 100000 minitagged_sorted.bam
```

## 3. Cell allele counting 
Requires vartrix https://github.com/10XGenomics/vartrix
and add /path/to/vartrix to your PATH
```
vartrix --umi --mapq 30 -b <bam file> -c <barcode tsv> --scoring-method coverage --threads 8 --ref-matrix ref.mtx --out-matrix alt.mtx -v <freebayes vcf> --fasta <fasta file used for remapping>
```
note the --threads argument and use an appropriate number of threads for your system.

## 4. Clustering cells by genotype
Requires Python3 with modules argparse, numpy, tensorflow
tensorflow requires Glibc >= 2.14
```
pip install tensorflow
```
should work if the glibc is up to date.
And go ahead and put the souporcell direcotry on your path
```
./souporcell.py -h
usage: cellection.py [-h] -a ALT_MATRIX -r REF_MATRIX -k NUM_CLUSTERS
                     [-l MAX_LOCI] [--min_alt MIN_ALT] [--min_ref MIN_REF]
                     [-t THREADS]

single cell RNAseq mixed genotype clustering using sparse mixture model
clustering with tensorflow.

optional arguments:
  -h, --help            show this help message and exit
  -a ALT_MATRIX, --alt_matrix ALT_MATRIX
                        alt matrix output from vartrix in coverage mode
  -r REF_MATRIX, --ref_matrix REF_MATRIX
                        ref matrix output from vartrix in coverage mode
  -k NUM_CLUSTERS, --num_clusters NUM_CLUSTERS
                        number of clusters to generate
  -l MAX_LOCI, --max_loci MAX_LOCI
                        maximum loci to consider per cell
  --min_alt MIN_ALT     minimum number of cells expressing the alt allele to
                        use the locus for clustering
  --min_ref MIN_REF     minimum number of cells expressing the ref allele to
                        use the locus for clustering
  -t THREADS, --threads THREADS
                        number of threads to run on
```
So generally something along the lines of
```
./souporcell.py -a alt.mtx -r ref.mtx -k <num_clusters> -t 8
```

## 5. Calling doublets
Rust required. To install rust:
```
curl https://sh.rustup.rs -sSf | sh
```
And to build troublet:
```
cd /path/to/souporcell/troublet
cargo build --release
```
And add /path/to/souporcell/troublet/target/release to your path
The usage is
```
/troublet -h
troublet 1.0
Haynes Heaton <whheaton@gmail.com>
Intergenotypic doublet detection given cluster assignments and cell allele counts

USAGE:
    troublet [OPTIONS] --alts <alts> --clusters <clusters>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -a, --alts <alts>                      alt allele counts per cell in sparse matrix format out of vartrix
    -c, --clusters <clusters>              cluster file output from schism
    -b, --debug <debug>...                 print debug info for index of cells listed
    -d, --doublet_prior <doublet_prior>    prior on doublets. Defaults to 0.5
    -r, --refs <refs>                      ref allele counts per cell in sparse matrix format out of vartrix
```
So generally
```
troublet -a alt.mtx -r ref.mtx --clusters clusters.tsv > doublet_output.tsv
```

## 6. Genotype and ambient RNA coinference
Python3 required with modules pystan, pyvcf, pickle, math, scipy, gzip (pip install should work for each)
```
consensus.py -h
usage: consensus.py [-h] -c CLUSTERS -a ALT_MATRIX -r REF_MATRIX [-p PLOIDY]
                    -d DOUBLETS -v VCF

consensus genotype calling and ambient RNA estimation

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        tsv cluster file
  -a ALT_MATRIX, --alt_matrix ALT_MATRIX
                        alt matrix file
  -r REF_MATRIX, --ref_matrix REF_MATRIX
                        ref matrix file
  -p PLOIDY, --ploidy PLOIDY
                        ploidy, must be 1 or 2, defaults to 2
  -d DOUBLETS, --doublets DOUBLETS
                        doublet calls
  -v VCF, --vcf VCF     vcf file from which alt and ref matrix were created
```
So generally
```
consensus.py -c clusters.tsv -a alt.mtx -r ref.mtx -d doublet_ouput.tsv --vcf <freebayes vcf>
```



