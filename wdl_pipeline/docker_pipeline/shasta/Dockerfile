FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

RUN apt-get install -y ncbi-blast+ python3 python3-pip && \
    pip3 install pybind11

WORKDIR /opt/shasta/
RUN wget https://github.com/chanzuckerberg/shasta/archive/0.3.0.tar.gz && \
    tar xvf 0.3.0.tar.gz && \
    cd shasta-0.3.0 && \
    mkdir build && \
    cd build && \
    ../scripts/InstallPrerequisites-Ubuntu.sh && \
    cmake .. && \
    make all -j 8 && \
    make install/strip && \
    mkdir /root/bin/shasta_0.3.0 && \
    cp /opt/shasta/shasta-0.3.0/build/shasta-install/bin/* /root/bin/shasta_0.3.0 && \
    mkdir /root/modules/shasta && \
    echo "#%Module" >>/root/modules/shasta/0.3.0 && \
    echo "append-path PATH /root/bin/shasta_0.3.0" >>/root/modules/shasta/0.3.0 && \
    echo "#%Module" >>/root/modules/shasta/.modulerc && \
    echo "module-version /0.3.0 default" >>/root/modules/shasta/.modulerc && \
    bash -i -c "module load shasta && shasta -h"

WORKDIR /opt/shasta/
RUN wget https://github.com/chanzuckerberg/shasta/archive/bf757a319bdd1dcca38ed86ae4d1dcf0f2a23f63.tar.gz && \
    tar xvf bf757a319bdd1dcca38ed86ae4d1dcf0f2a23f63.tar.gz && \
    cd shasta-bf757a319bdd1dcca38ed86ae4d1dcf0f2a23f63 && \
    mkdir build && \
    cd build && \
    ../scripts/InstallPrerequisites-Ubuntu.sh && \
    cmake .. && \
    make all -j 8 && \
    make install/strip && \
    mkdir /root/bin/shasta_0.3.0_bf757a3 && \
    cp /opt/shasta/shasta-bf757a319bdd1dcca38ed86ae4d1dcf0f2a23f63/build/shasta-install/bin/* /root/bin/shasta_0.3.0_bf757a3 && \
    echo "#%Module" >>/root/modules/shasta/0.3.0_bf757a3 && \
    echo "append-path PATH /root/bin/shasta_0.3.0_bf757a3" >>/root/modules/shasta/0.3.0_bf757a3 && \
    bash -i -c "module load shasta/0.3.0_bf757a3 && shasta -h"

#WORKDIR /root/scripts/shasta/
#COPY tmp/*.sh ./
WORKDIR /data
