#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


WORKDIR /root/scripts/bionano/
COPY _docker_bionano.sh tmp/* ./
COPY trf409.linux64 tmp/
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/bionano/_docker_bionano.sh"]