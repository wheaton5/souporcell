Bootstrap: docker

From: conda/miniconda3

%environment
	PATH=/opt/bedtools2/bin:/usr/local/envs/py36/bin:/opt/souporcell:/opt/souporcell/troublet/target/release:/usr/local/condabin:/opt/minimap2-2.7:/root/.cargo/bin:/opt/vartrix-v1.1.3-x86_64-linux/:/opt:$PATH

%post -c /bin/bash
        apt update
        yes | apt-get install wget
        yes | apt-get install build-essential
	yes | apt-get install curl
	echo blah
	yes | /usr/local/bin/conda create -n py36 python=3.6
	. /opt/conda/bin/activate py36
        yes | apt-get install libncurses5-dev
        yes | apt-get install zlib1g-dev
        yes | apt-get install libbz2-dev
        yes | apt-get install liblzma-dev
	cd /opt
	wget https://github.com/lh3/minimap2/archive/v2.7.tar.gz
        tar -xzvf v2.7.tar.gz
        cd minimap2-2.7
        make
	cd ..
	wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
	tar -zxvf bedtools-2.28.0.tar.gz
	cd bedtools2
	make
	cd ..
	CARGO_HOME=/opt/.cargo RUSTUP_HOME=/opt/.cargo bash -c 'curl https://sh.rustup.rs -sSf | sh -s -- -y'
	source /opt/.cargo/env
	which cargo
	rustup default stable
	yes | apt-get install git
	cd /opt
	git clone https://github.com/wheaton5/souporcell.git
	cd souporcell/troublet
	cargo build --release
	cd /opt/souporcell/souporcell
	cargo build --release
	cd /opt
	yes | /usr/local/envs/py36/bin/pip install pysam
	/usr/local/envs/py36/bin/pip install pyvcf
        /usr/local/envs/py36/bin/pip install numpy
        /usr/local/envs/py36/bin/pip install scipy
        /usr/local/envs/py36/bin/pip install pystan==2.17.1.0
        /usr/local/envs/py36/bin/pip install pyfaidx
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
	wget https://github.com/ekg/freebayes/releases/download/v1.3.1/freebayes-v1.3.1
	mv freebayes-v1.3.1 freebayes
	chmod 777 freebayes
        wget https://github.com/10XGenomics/vartrix/releases/download/v1.1.16/vartrix_linux
        mv vartrix_linux vartrix
        chmod 777 vartrix
