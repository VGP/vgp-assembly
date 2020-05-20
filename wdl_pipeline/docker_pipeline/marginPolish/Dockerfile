FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

WORKDIR /opt/marginPolish/
RUN git clone https://github.com/UCSC-nanopore-cgl/MarginPolish.git MarginPolish_1.2.0 && \
    cd MarginPolish_1.2.0 && \
    git fetch && \
    git checkout v1.2.0 && \
    git submodule update --init && \
    rm -rf /opt/marginPolish/MarginPolish_1.2.0/.git && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 8 && \
    mkdir /root/bin/marginPolish_1.2.0 && \
    cp /opt/marginPolish/MarginPolish_1.2.0/build/marginPolish /root/bin/marginPolish_1.2.0 && \
    mkdir /root/modules/marginPolish && \
    echo "#%Module" >>/root/modules/marginPolish/1.2.0 && \
    echo "append-path PATH /root/bin/marginPolish_1.2.0" >>/root/modules/marginPolish/1.2.0 && \
    echo "#%Module" >>/root/modules/marginPolish/.modulerc && \
    echo "module-version /1.2.0 default" >>/root/modules/marginPolish/.modulerc && \
    bash -i -c "module load marginPolish && marginPolish -h"

#WORKDIR /root/scripts/shasta/
#COPY tmp/*.sh ./
WORKDIR /data
