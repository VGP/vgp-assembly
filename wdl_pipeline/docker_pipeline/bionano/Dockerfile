FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

#TODO actually install bionano

WORKDIR /root/scripts/bionano/
COPY tmp/*.sh ./

WORKDIR /root/config/bionano
COPY tmp/*.xml ./

WORKDIR /data
