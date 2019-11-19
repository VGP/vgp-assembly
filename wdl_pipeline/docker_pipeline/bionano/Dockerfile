FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

WORKDIR /root/scripts/bionano/
COPY _docker_bionano.sh tmp/* ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/bionano/_docker_bionano.sh"]
