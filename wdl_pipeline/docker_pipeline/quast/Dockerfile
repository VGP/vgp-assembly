FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

#RUN apt-get update && apt-get install -y zlib1g-dev pkg-config libfreetype6-dev libpng-dev python-matplotlib python-setuptools libsasl2-dev libssl-dev

### QUAST
WORKDIR /root/tools/quast
RUN wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz && \
    tar -xzf quast-5.0.2.tar.gz && \
    rm quast-5.0.2.tar.gz && \
    cd quast-5.0.2 && \
    bash -i -c "module load python/2.7 && \
    which python && \
    python setup.py install"

WORKDIR /data
