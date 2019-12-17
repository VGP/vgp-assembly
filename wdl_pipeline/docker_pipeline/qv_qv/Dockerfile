FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


WORKDIR /root/scripts/qv/
COPY _docker_qv.sh tmp/*.sh tmp/*.jar ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/qv/_docker_qv.sh"]
