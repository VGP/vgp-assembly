#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# install cmake
WORKDIR /opt/cmake_install
RUN mkdir /opt/cmake && \
    wget https://cmake.org/files/v3.13/cmake-3.13.5-Linux-x86_64.sh && \
    sh cmake-3.13.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake && \
    rm -r /opt/cmake_install

WORKDIR /opt/bamtools
RUN wget https://github.com/pezmaster31/bamtools/archive/v2.5.1.tar.gz && \
    tar -xvf v2.5.1.tar.gz && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/bamtools bamtools-2.5.1/ && \
    make install && \
    mkdir /root/bin/bamtools_2.5.1 && \
    cp bin/bamtools /root/bin/bamtools_2.5.1 && \
    rm -r /opt/bamtools && \
    mkdir /root/modules/bamtools && \
    echo "#%Module" >>/root/modules/bamtools/2.5.1 && \
    echo "append-path PATH /root/bin/bamtools_2.5.1" >>/root/modules/bamtools/2.5.1 && \
    echo "#%Module" >>/root/modules/bamtools/.modulerc && \
    echo "module-version /2.5.1 default" >>/root/modules/bamtools/.modulerc

WORKDIR /opt/blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    mkdir /root/bin/blast_2.2.30+ && \
    mv /opt/blast/ncbi-blast-2.2.30+/bin/* /root/bin/blast_2.2.30+/ && \
    rm -r /opt/blast && \
    mkdir /root/modules/blast && \
    echo "#%Module" >>/root/modules/blast/2.2.30+ && \
    echo "append-path PATH /root/bin/blast_2.2.30+" >>/root/modules/blast/2.2.30+ && \
    echo "#%Module" >>/root/modules/blast/.modulerc && \
    echo "module-version /2.2.30+ default" >>/root/modules/blast/.modulerc

WORKDIR /opt/hmmer
RUN wget https://github.com/EddyRivasLab/hmmer/archive/hmmer-3.2.1.tar.gz && \
    tar -xvf hmmer-3.2.1.tar.gz && \
    cd hmmer-hmmer-3.2.1 && \
    git clone https://github.com/EddyRivasLab/easel.git && \
    cd easel && \
    git checkout 924d7efff9a765d1d5807904f9686e2513a1e4f2 && \
    cd .. && \
    autoconf && \
    ./configure --prefix /opt/hmmer/hmmer-hmmer-3.2.1/install/ && \
    make && \
    make install && \
    mkdir /root/bin/hmmer_3.2.1 && \
    cp install/bin/* /root/bin/hmmer_3.2.1 && \
    rm -r /opt/hmmer && \
    mkdir /root/modules/hmmer && \
    echo "#%Module" >>/root/modules/hmmer/3.2.1 && \
    echo "append-path PATH /root/bin/hmmer_3.2.1" >>/root/modules/hmmer/3.2.1 && \
    echo "#%Module" >>/root/modules/hmmer/.modulerc && \
    echo "module-version /3.2.1 default" >>/root/modules/hmmer/.modulerc

WORKDIR /opt/R
RUN apt-get install -y gfortran build-essential fort77 xorg-dev libblas-dev gcc-multilib gobjc++ aptitude libreadline-dev
RUN apt-get install -y libpcre3 libpcre3-dev
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz && \
    tar xvf R-3.2.1.tar.gz && \
    cd R-3.2.1 && \
    ./configure --enable-utf8 && \
    make && \
    mkdir /root/bin/R_3.2.1 && \
    cp -r bin/* /root/bin/R_3.2.1 && \
    rm -r /opt/R && \
    mkdir /root/modules/R && \
    echo "#%Module" >>/root/modules/R/3.2.1 && \
    echo "append-path PATH /root/bin/R_3.2.1" >>/root/modules/R/3.2.1 && \
    echo "#%Module" >>/root/modules/R/.modulerc && \
    echo "module-version /3.2.1 default" >>/root/modules/R/.modulerc

WORKDIR /opt/python3
RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
    tar xvf Python-3.6.0.tgz && \
    cd Python-3.6.0 && \
    ./configure && \
    make && \
    mkdir /root/bin/python_3.6.0 && \
    cp python /root/bin/python_3.6.0 && \
    rm -r /opt/python3 && \
    mkdir /root/modules/python && \
    echo "#%Module" >>/root/modules/python/3.6 && \
    echo "append-path PATH /root/bin/python_3.6.0" >>/root/modules/python/3.6 && \
    echo "#%Module" >>/root/modules/python/.modulerc && \
    echo "module-version /3.6 default" >>/root/modules/python/.modulerc



# TODO
#module load python/3.6

WORKDIR /root/scripts/busco/
COPY _docker_busco.sh tmp/*.sh ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/busco/_docker_busco.sh"]
