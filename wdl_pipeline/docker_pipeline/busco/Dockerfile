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

# TODO
#module load blast/2.2.30+	
#module load hmmer
#module load R
#module load python/3.6

WORKDIR /root/scripts/busco/
COPY _docker_busco.sh tmp/*.sh ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/busco/_docker_busco.sh"]
