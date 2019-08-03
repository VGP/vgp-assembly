# Makefile for bionano enzyme 2 asset - edited from similar file authored by Yih-Chii Hwang

SHELL=/bin/bash -e -x -o pipefail

UBUNTU_VERSION := $(shell test -e /etc/issue && (grep -o "Ubuntu [0-9]\+\.[0-9]\+" /etc/issue | sed -e "s/Ubuntu //"))

ifeq ($(UBUNTU_VERSION), 12.04)
        CRAN_SITE = "deb http://cran.rstudio.com/bin/linux/ubuntu precise/"
else ifeq ($(UBUNTU_VERSION), 14.04)
       CRAN_SITE="deb http://cran.rstudio.com/bin/linux/ubuntu trusty/"
endif

SRC="/etc/apt/sources.list"

IS_FOUND = $(shell sudo grep $(CRAN_SITE) $(SRC))
TMP_DIR="/tmp"
ZLIB_URL="http://zlib.net/zlib-1.2.11.tar.gz"
ZLIB_FILE="zlib-1.2.11.tar.gz"
ZLIB_DIR="zlib-1.2.11"

CURL_URL="https://curl.haxx.se/download/curl-7.54.1.tar.gz"
CURL_FILE="curl-7.54.1.tar.gz"
CURL_DIR="curl-7.54.1"

XORG_URL="http://archive.ubuntu.com/ubuntu/pool/main/x/xorg/xorg_7.7+1ubuntu8.1.tar.gz"
XORG_FILE="xorg_7.7+1ubuntu8.1.tar.gz"
XORG_DIR="xorg_7.7+1ubuntu8.1"

R_URL="https://cran.r-project.org/src/base/R-3/R-3.5.0.tar.gz"
R_FILE="R-3.5.0.tar.gz"
R_DIR="R-3.5.0"

LIBXML_URL="ftp://xmlsoft.org/libxml2/libxml2-2.7.2.tar.gz"
LIBXML_FILE="libxml2-2.7.2.tar.gz"
LIBXML_DIR="libxml2-2.7.2"

# bundle url created from file-FFVZ0KPKGbk5pyPxJGpvQ1FJ
BIONANO_BUNDLE_URL="https://dl.dnanex.us/F/D/7FzGkPzbK22339V0yP2V5j7g24xbXj5q1JFVk231/Solve3.2.1_04122018.tar.gz"
BIONANO_BUNDLE="Solve3.2.1_04122018.tar.gz"
SCAFFOLD_DIR="/Solve3.2.1_04122018/HybridScaffold/04122018"

all:

ifeq ($(IS_FOUND), )
	$(shell sudo echo $(CRAN_SITE) >> $(SRC))
endif
	# Trust the signing key for this repo.Reference: http://cran.rstudio.com/bin/linux/ubuntu/README.html
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	#Update apt-get
	sudo apt-get --yes update
	sudo apt-get --yes upgrade
	sudo apt-get --yes install xorg-dev 
	sudo apt-get --yes install zlib1g-dev liblzma-dev
	
	sudo curl -L $(R_URL) > $(TMP_DIR)/$(R_FILE)
	sudo tar -zxf $(TMP_DIR)/$(R_FILE) -C $(TMP_DIR)
	
	echo "Downloading and installing zlib-1.2.11........."
	sudo curl -L $(ZLIB_URL) > $(TMP_DIR)/$(ZLIB_FILE)
	sudo tar -zxf $(TMP_DIR)/$(ZLIB_FILE) -C $(TMP_DIR)
	cd $(TMP_DIR)/$(ZLIB_DIR); ./configure
	sudo $(MAKE) -C $(TMP_DIR)/$(ZLIB_DIR)
	sudo $(MAKE) -C $(TMP_DIR)/$(ZLIB_DIR) install

	echo "Downloading and installing curl-7.54.1........."
	sudo curl -L $(CURL_URL) > $(TMP_DIR)/$(CURL_FILE)
	sudo tar -zxf $(TMP_DIR)/$(CURL_FILE) -C $(TMP_DIR)
	cd $(TMP_DIR)/$(CURL_DIR); ./configure --with-ssl
	sudo $(MAKE) -C $(TMP_DIR)/$(CURL_DIR)
	sudo $(MAKE) -C $(TMP_DIR)/$(CURL_DIR) install

	echo "Downloading and installing R-3.5.0........."
	cd $(TMP_DIR)/$(R_DIR); ./configure
	sudo $(MAKE) -C $(TMP_DIR)/$(R_DIR)
	sudo $(MAKE) -C $(TMP_DIR)/$(R_DIR) install

	echo "Downloading libxml...."
	sudo apt-get install libxml2-dev

	echo "Installing all R dependencies........"
	sudo R --quiet -e "install.packages('data.table', repos='http://cran.rstudio.com/'); install.packages('igraph',repos='http://cran.rstudio.com/'); install.packages('intervals', repos='http://cran.rstudio.com/'); install.packages('MASS',repos='http://cran.rstudio.com/'); install.packages('parallel',repos='http://cran.rstudio.com/'); install.packages('XML',repos='http://cran.rstudio.com/'); install.packages('argparser',repos='http://cran.rstudio.com/'); install.packages('httpuv',repos='http://cran.rstudio.com/'); install.packages('shiny',repos='http://cran.rstudio.com/');"
	
	echo "Installing Bionano tools......"
	wget $(BIONANO_BUNDLE_URL) -q
	sudo tar zxvf $(BIONANO_BUNDLE) -C /
	Rscript $(SCAFFOLD_DIR)/runTGH.R --help

