FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

### Salsa
# 3.0.2
WORKDIR /root/tools/salsa
RUN wget https://github.com/machinegun/SALSA/archive/v2.2.tar.gz && \
    tar -zxf v2.2.tar.gz && \
    rm v2.2.tar.gz && \
    cd SALSA-2.2 && \
    make

WORKDIR /root/scripts/salsa/
COPY tmp/*.sh ./
WORKDIR /data
