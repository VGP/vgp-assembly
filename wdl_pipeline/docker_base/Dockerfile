FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ARG git_commit

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ vim tcl8.6-dev autoconf autopoint python-dev \
                       perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# get modules
WORKDIR /opt
RUN git clone https://github.com/cea-hpc/modules.git && \
    cd /opt/modules && \
    git fetch && \
    git checkout $git_commit 

# install modules
WORKDIR /opt/modules
RUN ./configure --with-tclx-ver=8.6 --with-tcl-ver=8.6 && \
    make && \
    make install

# user configuration
RUN echo "source /usr/local/Modules/init/bash" >>/root/.bashrc && \
    echo "module use /root/modules/" >>/root/.bashrc && \
    mkdir /root/modules && \
    mkdir /root/bin

### minimap2
# 2.11 
WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/archive/v2.11.tar.gz && \
    tar -xvf v2.11.tar.gz && \
    cd minimap2-2.11 && \
    make && \
    mkdir /root/bin/minimap2_2.11 && \
    cp minimap2 /root/bin/minimap2_2.11 && \
    cd /opt/minimap2 && \
    rm -r /opt/minimap2/* && \
    mkdir /root/modules/minimap2 && \
    echo "#%Module" >>/root/modules/minimap2/2.11 && \
    echo "append-path PATH /root/bin/minimap2_2.11" >>/root/modules/minimap2/2.11 && \
    echo "#%Module" >>/root/modules/minimap2/.modulerc && \
    echo "module-version /2.11 default" >>/root/modules/minimap2/.modulerc

### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
    cd samtools-1.9/ && \ 
    autoheader && \
    autoconf -Wno-header && \
    ./configure --without-curses && \
    make && \
    mkdir /root/bin/samtools_1.9 && \
    cp samtools /root/bin/samtools_1.9 && \
    rm -r /opt/samtools && \
    mkdir /root/modules/samtools && \
    echo "#%Module" >>/root/modules/samtools/1.9 && \
    echo "append-path PATH /root/bin/samtools_1.9" >>/root/modules/samtools/1.9 && \
    echo "#%Module" >>/root/modules/samtools/.modulerc && \
    echo "module-version /1.9 default" >>/root/modules/samtools/.modulerc 

