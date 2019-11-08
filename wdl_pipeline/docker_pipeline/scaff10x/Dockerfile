#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


WORKDIR /root/tools/scaff10x/
RUN git clone https://github.com/wtsi-hpag/Scaff10X.git Scaff10X-4.1 && \
    cd Scaff10X-4.1 && \
    git fetch && \
    git checkout v4.1 && \
    rm -rf .git && \
    ./install.sh

WORKDIR /root/scripts/scaff10x/
COPY tmp/*.sh ./
WORKDIR /data
