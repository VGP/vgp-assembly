#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


WORKDIR /root/scripts/qv/
COPY _docker_qv.sh tmp/*.sh tmp/*.jar ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/qv/_docker_qv.sh"]
