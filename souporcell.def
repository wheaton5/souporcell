BootStrap: library
From: ubuntu:18.04

%environment
    # set up all essential environment variables
    export LC_ALL=C
    export PATH=/miniconda3/bin:/opt/minimap2-2.26_x64-linux:/opt/bedtools-2:/opt/freebayes/build:/opt/hisat2-2.2.1:/opt/bedtools2/bin:/opt/souporcell:/opt/souporcell/troublet/target/release:/usr/local/condabin:/opt/minimap2-2.26_x64-linux:/root/.cargo/bin:/opt:$PATH
    export PYTHONPATH=/miniconda3/lib/python3.9/:$PYTHONPATH
    
    # activate conda environment
    source activate base;
    conda activate;
    
%post
    # update and install essential dependencies
    apt-get -y update
    apt-get update && apt-get install -y automake build-essential bzip2 wget git default-jre unzip
    apt-get -y install wget
    apt-get -y install curl
    
    # download, install, and update miniconda3
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /miniconda3/
    rm Miniconda3-latest-Linux-x86_64.sh
    
    # install dependencies via conda
    export PATH="/miniconda3/bin:$PATH"
    conda install -y -c conda-forge pip numpy # general dependencies
    conda update --all
    
    apt-get -y install libncurses5-dev
    apt-get -y install zlib1g-dev
    apt-get -y install libbz2-dev
    apt-get -y install liblzma-dev
    apt-get -y install zip unzip
    apt-get -y install pkg-config
    cd /opt
    wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
    mv download hisat2.zip
    unzip hisat2.zip
    wget https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2
    tar -xf minimap2-2.26_x64-linux.tar.bz2
    cd ..
    wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
    tar -xf bedtools-2.28.0.tar.gz
    cd bedtools2
    make
    cd ..
    CARGO_HOME=/opt/.cargo RUSTUP_HOME=/opt/.cargo bash -c 'curl https://sh.rustup.rs -sSf | sh -s -- -y'
    . /opt/.cargo/env
    which cargo
    rustup default stable
    apt-get -y install git
    cd /opt
    git clone --recursive https://github.com/wheaton5/souporcell.git
    cd souporcell/troublet
    cargo build --release
    cd /opt/souporcell/souporcell
    cargo build --release
    cd /opt
    pip install pysam
    pip install git+https://github.com/stan-dev/pystan2.git@master
    pip install pyfaidx
    pip install "setuptools<58" --upgrade
    pip install pyvcf
    pip install scipy
    cd /opt
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
    tar xvfj htslib-1.9.tar.bz2
    cd htslib-1.9
    ./configure
    make
    make install
    cd ..
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar xvfj samtools-1.9.tar.bz2
    rm samtools-1.9.tar.bz2
    cd samtools-1.9
    ./configure
    make
    make install
    cd ..
    wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
    tar xvfj bcftools-1.9.tar.bz2
    rm bcftools-1.9.tar.bz2
    cd bcftools-1.9
    ./configure
    make
    make install
    cd /opt
    wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
    gunzip freebayes-1.3.6-linux-amd64-static.gz
    mv freebayes-1.3.6-linux-amd64-static freebayes
    chmod 777 freebayes
    wget https://github.com/10XGenomics/vartrix/releases/download/v1.1.22/vartrix_linux
    mv vartrix_linux vartrix
    chmod 777 vartrix
    cd ..
    

%labels
    Author whheaton
    Version v2.1
    MyLabel souporcell

