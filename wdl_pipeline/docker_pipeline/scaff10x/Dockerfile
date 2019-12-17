FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


WORKDIR /root/tools/scaff10x/
RUN git clone https://github.com/wtsi-hpag/Scaff10X.git Scaff10X-4.1 && \
    cd Scaff10X-4.1 && \
    git fetch && \
    git checkout v4.1 && \
    rm -rf .git && \
    ./install.sh

RUN apt-get install -y libgmp-dev libmpfr-dev libmpc-dev
WORKDIR /opt/gcc
RUN wget http://www.netgull.com/gcc/releases/gcc-7.4.0/gcc-7.4.0.tar.gz && \
    tar xvf gcc-7.4.0.tar.gz && \
    cd gcc-7.4.0 && \
    ./contrib/download_prerequisites && \
    cd .. && \
    mkdir build_7.4.0 && \
    cd build_7.4.0 && \
    ../gcc-7.4.0/configure --prefix=/opt/gcc/build_7.4.0 --enable-languages=c,c++,fortran && \
    make -j 12 && \
    make install && \
    mkdir /root/bin/gcc_7.4.0  && \
	cp -r /opt/gcc/build_7.4.0/bin/* /root/bin/gcc_7.4.0 && \
	rm -r /opt/gcc && \
	mkdir /root/modules/gcc && \
	echo "#%Module" >>/root/modules/gcc/7.4.0 && \
    echo "prepend-path PATH /root/bin/gcc_7.4.0" >>/root/modules/gcc/7.4.0 && \
    echo "#%Module" >>/root/modules/gcc/.modulerc && \
    echo "module-version /7.4.0 default" >>/root/modules/gcc/.modulerc

WORKDIR /root/scripts/scaff10x/
COPY tmp/*.sh ./
WORKDIR /data
