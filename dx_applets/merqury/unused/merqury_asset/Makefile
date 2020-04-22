SHELL=/bin/bash -e -x -o pipefail
TMP_DIR="/tmp"
R_URL="https://cran.r-project.org/src/base/R-3/R-3.5.0.tar.gz"
R_FILE="R-3.5.0.tar.gz"
R_DIR="R-3.5.0"
CRAN_SITE="deb http://cran.rstudio.com/bin/linux/ubuntu trusty/"
all:
	sudo apt-get -y update && apt-get -y upgrade 	
	# Trust the signing key for this repo.Reference: http://cran.rstudio.com/bin/linux/ubuntu/README.html
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# install xorg-dev for R
	sudo apt-get -y install xorg-dev 
	# git parallel and common bioinfo tools
	sudo apt-get -y install git bedtools ncbi-blast+ parallel
	# for ubuntu14.04
	sudo apt-get -y install liblzma-dev
	# HTSlib
	sudo wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar -xjf htslib-1.9.tar.bz2
	pushd htslib-1.9 && sudo ./configure && sudo make && sudo make install && popd
	# samtools	
	sudo wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar xfj samtools-1.9.tar.bz2
	pushd samtools-1.9 && sudo ./configure && sudo make && sudo make install && popd
	# Rscript
	echo "Downloading and installing R-3.5.0........."
	sudo curl -L $(R_URL) > $(TMP_DIR)/$(R_FILE)
	sudo tar -zxf $(TMP_DIR)/$(R_FILE) -C $(TMP_DIR)
	cd $(TMP_DIR)/$(R_DIR); ./configure
	sudo $(MAKE) -C $(TMP_DIR)/$(R_DIR)
	sudo $(MAKE) -C $(TMP_DIR)/$(R_DIR) install	
	sudo su - -c "R -e \"install.packages('ggplot2', repos='http://cran.rstudio.com/')\""
	sudo su - -c "R -e \"install.packages('scales', repos='http://cran.rstudio.com/')\""
	# lastz	
	export LASTZ_INSTALL=/usr/bin
	sudo wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.00.tar.gz
	sudo tar xzpf lastz-1.04.00.tar.gz
	pushd lastz-distrib-1.04.00/ && sudo make && sudo make install && popd
	sudo mv lastz-distrib-1.04.00/bin/* /usr/bin
	# minimap2
	sudo wget https://github.com/lh3/minimap2/releases/download/v2.12/minimap2-2.12.tar.bz2 && tar -xjf minimap2-2.12.tar.bz2
	pushd minimap2-2.12 && sudo make && sudo mv minimap2 /usr/bin && popd
	# mummer4
	sudo wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
	tar -zxvf mummer-4.0.0beta2.tar.gz
	pushd mummer-4.0.0beta2 && ./configure --enable-python-binding && sudo make && sudo make install && popd; 
	# meryl
	meryl-1.0.Linux-amd64.tar.xz

