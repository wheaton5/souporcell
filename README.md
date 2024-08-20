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
```
singularity pull --arch amd64 library://wheaton5/souporcell/souporcell:release
```

If you are running on a scientific cluster, they will likely have singularity, contact your sysadmin for more details. 
If you are running on your own linux box you may need to install [singularity](https://www.sylabs.io/guides/3.2/user-guide/quick_start.html#quick-installation-steps)

requires singularity >= 3.0 or apptainer version >= 1.0.0 (singularity rebranded as apptainer and changed its version numbers)
```
which singularity
singularity --version
```
You should now be able to run souporcell_pipeline.py through the singularity container. Singularity automatically mounts the current working directory and directories downstream from where you run it, otherwise you would need to manually mount those directories. Therefor you can run it from a directory that is upstream of all of the inputs. Input files are the cellranger bam, cellranger barcodes file, and a reference fasta. The cellranger bam is located in the cellranger outs directory and is called possorted_genome_bam.bam. The barcodes file is located in the cellranger outs/filtered_gene_bc_matrices/<ref_name>/barcodes.tsv. The reference fasta should be of the same species but does not necessarily need to be the exact cellranger reference. 

The options for using souporcell are:
```
singularity exec souporcell_latest.sif souporcell_pipeline.py -h
usage: souporcell_pipeline.py [-h] -i BAM -b BARCODES -f FASTA -t THREADS -o
                              OUT_DIR -k CLUSTERS [-p PLOIDY]
                              [--min_alt MIN_ALT] [--min_ref MIN_REF]
                              [--max_loci MAX_LOCI] [--restarts RESTARTS]
                              [--common_variants COMMON_VARIANTS]
                              [--known_genotypes KNOWN_GENOTYPES]
                              [--known_genotypes_sample_names KNOWN_GENOTYPES_SAMPLE_NAMES [KNOWN_GENOTYPES_SAMPLE_NAMES ...]]
                              [--skip_remap SKIP_REMAP] [--ignore IGNORE]

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
  --restarts RESTARTS   number of restarts in clustering, when there are > 12
                        clusters we recommend increasing this to avoid local
                        minima
                         --common_variants COMMON_VARIANTS
                        common variant loci or known variant loci vcf, must be
                        vs same reference fasta
  --known_genotypes KNOWN_GENOTYPES
                        known variants per clone in population vcf mode, must
                        be .vcf right now we dont accept gzip or bcf sorry
  --known_genotypes_sample_names KNOWN_GENOTYPES_SAMPLE_NAMES [KNOWN_GENOTYPES_SAMPLE_NAMES ...]
                        which samples in population vcf from known genotypes
                        option represent the donors in your sample
  --skip_remap SKIP_REMAP
                        don't remap with minimap2 (not recommended unless in
                        conjunction with --common_variants
  --ignore IGNORE       set to True to ignore data error assertions
```
A typical command looks like
```
singularity exec /path/to/souporcell_latest.sif souporcell_pipeline.py -i /path/to/possorted_genome_bam.bam -b /path/to/barcodes.tsv -f /path/to/reference.fasta -t num_threads_to_use -o output_dir_name -k num_clusters
```
The above command will run all six steps of the pipeline and it will require up to 24gb of ram for human (minimap2 bam index is high water mark for memory). For smaller genomes, fewer clusters, lower --max-loci will require less memory. Note that souporcell will require roughly 2x the amount of diskspace that the input bam file takes up. This dataset should take several hours to run on 8 threads mostly due to read processing, remapping, and variant calling.

If you have a common snps file you may want to use the --common_variants option with or without the --skip_remap option. This option will skip conversion to fastq, remapping with minimap2, and reattaching barcodes, and the --common_variants will remove the freebayes step. Each which will save a significant amount of time, but --skip-remap isn't recommended without --common_variants.

Common variant files from 1k genomes filtered to variants >= 2% allele frequency in the population and limited to SNPs can be found here for GRCh38
```
curl ftp://ftp.eng.auburn.edu/pub/whh0027/common_variants_grch38.vcf.gz -o common_variants_grch38.vcf.gz
```
or for hg19 here
```
curl ftp://ftp.eng.auburn.edu/pub/whh0027/filtered_2p_1kgenomes_hg19.vcf.gz -o common_variants_hg19.vcf.gz
```

## Practice/Testing data set: Demuxlet paper data
```
wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398235/A.merged.bam.1 -O A.merged.bam
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560245/suppl/GSM2560245_barcodes.tsv.gz
gunzip GSM2560245_barcodes.tsv.gz
```
And if you don't have a human reference sitting around, grab one here
```
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz
```
Now you should be ready to test it out
```
singularity exec /path/to/souporcell_latest.sif souporcell_pipeline.py -i A.merged.bam -b GSM2560245_barcodes.tsv -f refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -t 8 -o demux_data_test -k 4
```

This should require about 20gb of ram mostly because of the minimap2 indexing step. I might soon host an index and reference for human to make this less painful.

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
I suggest you use the conda env I have set up by using the following command if you have conda or miniconda
```
conda env create -f /path/to/souporcell/souporcell_env.yaml
conda activate souporcell
```
You will also need Rust and to compile the two rust binaries
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cd /path/to/souporcell/souporcell && cargo build --release
cd /path/to/souporcell/troublet && cargo build --release
```
otherwise python packages tensorflow, pyvcf, pystan, pyfaidx, numpy, scipy are required, but as the versions change, I do recommend using the presetup env.

## To run through the pipeline script
```
souporcell_pipeline.py -i /path/to/possorted_genome_bam.bam -b /path/to/barcodes.tsv -f /path/to/reference.fasta -t num_threads_to_use -o output_dir_name -k num_clusters
```

## To run things step by step not through the pipeline script

### 1. Remapping
We discuss the need for remapping in our manuscript. We need to keep track of cell barcodes and and UMIs, so we first create a fastq with those items encoded in the readname.
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
freebayes -f <reference_fasta> -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --limit-coverage 100000 minitagged_sorted.bam
```

### 3. Cell allele counting 
Requires [vartrix](https://github.com/10XGenomics/vartrix)
and add /path/to/vartrix to your PATH
```
vartrix --umi --mapq 30 -b <bam file> -c <barcode tsv> --scoring-method coverage --threads 8 --ref-matrix ref.mtx --out-matrix alt.mtx -v <freebayes vcf> --fasta <fasta file used for remapping>
```
note the --threads argument and use an appropriate number of threads for your system.

### 4. Clustering cells by genotype
Rust required. To install rust:
```
curl https://sh.rustup.rs -sSf | sh
```
and to build souporcell clustering
```
cd /path/to/souporcell/souporcell
cargo build --release
```
And add /path/to/souporcell/souporcell/target/release to your path
usage
```
souporcell -h
souporcell 2.4
Haynes Heaton <whheaton@gmail.com>
clustering scRNAseq cells by genotype

USAGE:
    souporcell [OPTIONS] --alt_matrix <alt_matrix> --barcodes <barcodes> --num_clusters <num_clusters> --ref_matrix <ref_matrix>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -a, --alt_matrix <alt_matrix>                                           alt matrix from vartrix
    -b, --barcodes <barcodes>                                               cell barcodes
        --initialization_strategy <initialization_strategy>
            cluster initialization strategy, defaults to kmeans++, valid values are kmeans++, random_uniform,
            middle_variance, random_cell_assignment
        --known_cell_assignments <known_cell_assignments>
            tsv with barcodes and their known assignments

    -g, --known_genotypes <known_genotypes>
            NOT YET IMPLEMENTED population vcf/bcf of known genotypes if available.
            
        --known_genotypes_sample_names <known_genotypes_sample_names>...
            NOT YET IMPLEMENTED sample names, must be samples from the known_genotypes vcf

        --min_alt <min_alt>
            minimum number of cells containing the alt allele for the variant to be used for clustering

        --min_alt_umis <min_alt_umis>                                       min alt umis to use locus for clustering
        --min_ref <min_ref>
            minimum number of cells containing the ref allele for the variant to be used for clustering

        --min_ref_umis <min_ref_umis>                                       min ref umis to use locus for clustering
    -k, --num_clusters <num_clusters>                                       number of clusters
    -r, --ref_matrix <ref_matrix>                                           ref matrix from vartrix
    -r, --restarts <restarts>                                               number of random seedings
        --seed <seed>                                                       optional random seed
    -t, --threads <threads>                                                 number of threads to use
```
So generally something along the lines of
```
souporcell -a alt.mtx -r ref.mtx -b barcodes.tsv -k <num_clusters> -t 8 > clusters_tmp.tsv
```
(note clusters_tmp.tsv output as the doublet caller outputs the final clusters file)

### 5. Calling doublets
Rust required.
Build troublet:
```
cd /path/to/souporcell/troublet
cargo build --release
```
And add /path/to/souporcell/troublet/target/release to your path
The usage is
```
troublet -h
troublet 2.4
Haynes Heaton <whheaton@gmail.com>
Intergenotypic doublet detection given cluster assignments and cell allele counts

USAGE:
    troublet [OPTIONS] --alts <alts> --clusters <clusters>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -a, --alts <alts>                              alt allele counts per cell in sparse matrix format out of vartrix
    -c, --clusters <clusters>                      cluster file output from schism
    -b, --debug <debug>...                         print debug info for index of cells listed
    -d, --doublet_prior <doublet_prior>            prior on doublets. Defaults to 0.5
        --doublet_threshold <doublet_threshold>    doublet posterior threshold, defaults to 0.90
    -r, --refs <refs>                              ref allele counts per cell in sparse matrix format out of vartrix
        --singlet_threshold <singlet_threshold>    singlet posterior threshold, defaults to 0.90
```
So generally
```
troublet -a alt.mtx -r ref.mtx --clusters clusters_tmp.tsv > clusters.tsv
```

### 6. Genotype and ambient RNA coinference
Python3 required with modules pystan, pyvcf, pickle, math, scipy, gzip (pip install should work for each)
```
consensus.py -h
usage: consensus.py [-h] -c CLUSTERS -a ALT_MATRIX -r REF_MATRIX [-p PLOIDY]
                    --soup_out SOUP_OUT --vcf_out VCF_OUT --output_dir
                    OUTPUT_DIR -v VCF

consensus genotype calling and ambient RNA estimation

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        tsv cluster file from the troublet output
  -a ALT_MATRIX, --alt_matrix ALT_MATRIX
                        alt matrix file
  -r REF_MATRIX, --ref_matrix REF_MATRIX
                        ref matrix file
  -p PLOIDY, --ploidy PLOIDY
                        ploidy, must be 1 or 2, defaults to 2
  --soup_out SOUP_OUT   soup output
  --vcf_out VCF_OUT     vcf output
  --output_dir OUTPUT_DIR
                        output directory
  -v VCF, --vcf VCF     vcf file from which alt and ref matrix were created
```
So generally
```
consensus.py -c clusters.tsv -a alt.mtx -r ref.mtx --soup_out soup.txt -v <freebayes vcf> --vcf_out cluster_genotypes.vcf --output_dir .
```



