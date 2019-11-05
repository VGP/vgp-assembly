#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV PURGE_DUPS_COMMIT=f58cd4dcc55d6aa5794ccd158c3be59b08430f27

### minimap2
# 2.11
WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/archive/v2.17.tar.gz && \
    tar -xvf v2.17.tar.gz && \
    rm -r /opt/minimap2/v2.17.tar.gz && \
    cd minimap2-2.17 && \
    make && \
    mkdir /root/bin/minimap2_2.17 && \
    cp minimap2 /root/bin/minimap2_2.17 && \
    echo "#%Module" >>/root/modules/minimap2/2.17 && \
    echo "append-path PATH /root/bin/minimap2_2.17" >>/root/modules/minimap2/2.17 && \
    echo "#%Module" >/root/modules/minimap2/.modulerc && \
    echo "module-version /2.17 default" >>/root/modules/minimap2/.modulerc

WORKDIR /root/tools/
RUN git clone https://github.com/dfguan/purge_dups.git && \
    cd purge_dups && \
    git fetch && \
    git checkout $PURGE_DUPS_COMMIT && \
    cd src && \
    make

WORKDIR /root/scripts/purge_dups/
COPY tmp/*.sh ./
WORKDIR /data
