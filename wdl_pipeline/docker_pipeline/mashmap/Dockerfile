FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


WORKDIR /root/scripts/mashmap/
COPY _docker_mashmap.sh tmp/*.sh ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/mashmap/_docker_mashmap.sh"]
