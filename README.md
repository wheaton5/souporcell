# souporcell 

<img src="https://github.com/wheaton5/souporcell/blob/master/souporcell_star.png" width="100">

Preprint manuscript of this method available at https://www.biorxiv.org/content/10.1101/699637v1

souporcell is a method for clustering mixed-genotype scRNAseq experiments by individual.

The inputs are just the possorted_genome_bam.bam, and barcodes.tsv as output from [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).
souporcell is comprised of 6 steps with the first 3 using external tools and the final using the code provided here.
1. Remapping ([minimap2](https://github.com/lh3/minimap2))
2. Calling candidate variants ([freebayes](https://github.com/ekg/freebayes))
3. Cell allele counting ([vartrix](https://github.com/10XGenomics/vartrix))
4. Clustering cells by genotype (souporcell.py)
5. Calling doublets (troublet)
6. Calling cluster genotypes and inferring amount of ambient RNA (consensus.py)

## Easy Installation (Linux) (recommended) 

Download singularity image (1.3gb) (singularity is similar to docker but safe for clusters)
Google drive makes it annoyingly difficult to download via the terminal. The following command will download and name souporcell.sif (singularity image file) to your current directory.
```
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1cbMuUJPwXd8BsUMY0XmtyWhu4Jj1ckOb' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1cbMuUJPwXd8BsUMY0XmtyWhu4Jj1ckOb" -O souporcell.sif && rm -rf /tmp/cookies.txt
```
```
ls -lah souporcell.sif
```
This file is ~1.3Gb, so please make sure you have enough space for it. You can also download manually [here](https://drive.google.com/open?id=1cbMuUJPwXd8BsUMY0XmtyWhu4Jj1ckOb)

If you are running on a scientific cluster, they will likely have singularity, contact your sysadmin for more details. 
If you are running on your own linux box you may need to install [singularity](https://www.sylabs.io/guides/3.2/user-guide/quick_start.html#quick-installation-steps)

requires singularity >= 3.0
```
which singularity
singularity --version
```
You should now be able to run souporcell_pipeline.py through the singularity container. Singularity automatically mounts the current working directory and directories downstream from where you run it, otherwise you would need to manually mount those directories. Therefor you can run it from a directory that is upstream of all of the inputs. Input files are the cellranger bam, cellranger barcodes file, and a reference fasta. The cellranger bam is located in the cellranger outs directory and is called possorted_genome_bam.bam. The barcodes file is located in the cellranger outs/filtered_gene_bc_matrices/<ref_name>/barcodes.tsv. The reference fasta should be of the same species but does not necessarily need to be the exact cellranger reference. 

The options for using souporcell are:
```
singularity exec souporcell.sif souporcell_pipeline.py -h
usage: souporcell_pipeline.py [-h] -i BAM -b BARCODES -f FASTA -t THREADS -o
                              OUT_DIR -k CLUSTERS [-p PLOIDY]
                              [--min_alt MIN_ALT] [--min_ref MIN_REF]
                              [--max_loci MAX_LOCI] [--ignore IGNORE]

single cell RNAseq mixed genotype clustering using sparse mixture model
clustering with tensorflow.

optional arguments:
  -h, --help            show this help message and exit
  -i BAM, --bam BAM     cellranger bam
  -b BARCODES, --barcodes BARCODES
                        barcodes.tsv from cellranger
  -f FASTA, --fasta FASTA
                        reference fasta file
  -t THREADS, --threads THREADS
                        max threads to use
  -o OUT_DIR, --out_dir OUT_DIR
                        name of directory to place souporcell files
  -k CLUSTERS, --clusters CLUSTERS
                        number cluster, tbd add easy way to run on a range of
                        k
  -p PLOIDY, --ploidy PLOIDY
                        ploidy, must be 1 or 2, default = 2
  --min_alt MIN_ALT     min alt to use locus, default = 10.
  --min_ref MIN_REF     min ref to use locus, default = 10.
  --max_loci MAX_LOCI   max loci per cell, affects speed, default = 2048.
  --ignore IGNORE       set to True to ignore data error assertions
```
A typical command looks like
```
singularity exec /path/to/souporcell.sif souporcell_pipeline.py -i /path/to/possorted_genome_bam.bam -b /path/to/barcodes.tsv -f /path/to/reference.fasta -t num_threads_to_use -o output_dir_name -k num_clusters
```
The recommended number of threads is 8. 

The above command will run all six steps of the pipeline and it will require up to 24gb of ram for human (minimap2 bam index is high water mark for memory). For smaller genomes, fewer clusters, lower --max-loci will require less memory. Note that souporcell will require roughly 2x the amount of diskspace that the input bam file takes up.

## Practice/Testing data set: Demuxlet paper data
```
wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398235/A.merged.bam.1 -O A.merged.bam
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560245/suppl/GSM2560245_barcodes.tsv.gz
gunzip GSM2560245_barcodes.tsv.gz
```
And if you don't have a human reference sitting around, grab one here
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
```
Now you should be ready to test it out
```
singularity exec /path/to/souporcell.sif souporcell_pipeline.py -i A.merged.bam -b GSM2560245_barcodes.tsv -f human_g1k_v37.fasta -t 8 -o demux_data_test -k 4
```

This should require about 20gb of ram mostly because of the minimap2 indexing step. I might soon host an index and reference for human to make this less painful.

Your output should look something like 
```
checking modules
imports done
checking fasta
creating chunks
generating fastqs with cell barcodes and umis in readname
remapping with minimap2
cleaning up tmp fastqs
repopulating cell barcode and UMI tags
cleaning up tmp samfiles
sorting retagged bam files
merging sorted bams
running freebayes
merging vcfs
running vartrix
running souporcell clustering
running souporcell doublet detection
running co inference of ambient RNA and cluster genotypes
/opt/conda/envs/py36/lib/python3.6/site-packages/pystan/misc.py:399: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  elif np.issubdtype(np.asarray(v).dtype, float):
Initial log joint probability = -56040.6
    Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes
       5      -29176.8    0.00143424     0.0384157           1           1       16
Optimization terminated normally:
  Convergence detected: relative gradient magnitude is below tolerance
```

and your output directory should have 
```
alt.mtx          cluster_genotypes.vcf  clusters.tsv  ref.mtx         souporcell_merged_sorted_vcf.vcf.gz      souporcell_minimap_tagged_sorted.bam
ambient_rna.txt  clusters_tmp.tsv       minimap.err   souporcell.log  souporcell_merged_sorted_vcf.vcf.gz.tbi  souporcell_minimap_tagged_sorted.bam.bai
```
The important files are 
1. clusters.tsv
2. cluster_genotypes.vcf
3. ambient_rna.txt

clusters.tsv will look like
```
barcode status  assignment      log_loss_singleton      log_loss_doublet        cluster0        cluster1
AAACCTGAGATCCGAG-1      singlet 0       -152.7778890920112      -190.5463095948822      -43.95302689281067      -101.63377524087669
AAACCTGAGCACCGTC-1      singlet 0       -78.56014177554212      -96.66255440088581      -20.949294849836267     -52.57478083591962
AAACCTGAGTACGATA-1      singlet 0       -216.0188863327174      -281.3888392065457      -63.059016939362536     -159.5450834682198
AAACCTGGTACATGTC-1      singlet 1       -47.189434469216565     -96.30865717225866      -62.652900832546955     -15.284168900754413
AAACCTGTCTACTCAT-1      singlet 0       -129.30104434183454     -167.87811467946756     -41.09158213888751      -106.3201962010145
AAACCTGTCTTGTCAT-1      singlet 0       -85.99781433701455      -110.81701038967158     -24.518165091815554     -60.05279033826837
AAACGGGCACTGTTAG-1      singlet 0       -154.26595878718032     -191.05465308213363     -31.356408693487197     -81.61186496254497
AAACGGGCATCATCCC-1      singlet 1       -46.33205678267174      -80.24152434540565      -50.78221280006256      -14.615983876840312
AAACGGGGTAGGGTAC-1      singlet 0       -240.5237900569412      -302.91575436035924     -71.79370547349878      -154.08594135029728
AAACGGGTCGGCATCG-1      singlet 0       -166.66827966974532     -226.56795157885028     -51.08790637893961      -148.04625123166286
```
With the cell barcode, singlet/doublet status, cluster, log_loss_singleton, log_loss_doublet, followed by log loss for each cluster.

2. cluster_genotypes.vcf is a vcf with genotypes for each cluster for each variant in the input vcf from freebayes

and 

3. ambient_rna.txt just contains the ambient RNA percentage detected

## Hard install

Instead of using singularity you can install everything independently (not recommended, but shouldn't be too bad)
```
git clone https://github.com/wheaton5/souporcell.git
```
put souporcell directory on your PATH 
requires samtools, bcftools, htslib, python3, freebayes, vartrix, minimap2 all on your PATH
python packages tensorflow, pyvcf, pystan, pyfasta, numpy, scipy

## To run through the pipeline script
```
souporcell_pipeline.py -i /path/to/possorted_genome_bam.bam -b /path/to/barcodes.tsv -f /path/to/reference.fasta -t num_threads_to_use -o output_dir_name -k num_clusters
```

## To run things step by step not through the pipeline script

### 1. Remapping
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
Requires [minimap2](https://github.com/lh3/minimap2)
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
Requires [samtools](http://www.htslib.org/)
```
samtools sort minitagged.bam minitagged_sorted.bam
samtools index minitagged_sorted.bam
```

### 2. Calling candidate variants
You may wish to break this into multiple jobs such as 1 job per chromosome and merge after but the basic command is the following.
Requires [freebayes](https://github.com/ekg/freebayes) and add /path/to/freebayes/bin to your PATH
```
freebayes -f <reference_fasta> -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --max-coverage 100000 minitagged_sorted.bam
```

### 3. Cell allele counting 
Requires [vartrix](https://github.com/10XGenomics/vartrix)
and add /path/to/vartrix to your PATH
```
vartrix --umi --mapq 30 -b <bam file> -c <barcode tsv> --scoring-method coverage --threads 8 --ref-matrix ref.mtx --out-matrix alt.mtx -v <freebayes vcf> --fasta <fasta file used for remapping>
```
note the --threads argument and use an appropriate number of threads for your system.

### 4. Clustering cells by genotype
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

### 5. Calling doublets
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

### 6. Genotype and ambient RNA coinference
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



