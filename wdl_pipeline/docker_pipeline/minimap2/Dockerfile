#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

WORKDIR /root/scripts/minimap2/
COPY _docker_minimap2.sh tmp/*.sh ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/minimap2/_docker_minimap2.sh"]
