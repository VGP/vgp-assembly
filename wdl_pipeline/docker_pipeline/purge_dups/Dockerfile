FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV PURGE_DUPS_COMMIT=f58cd4dcc55d6aa5794ccd158c3be59b08430f27

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
