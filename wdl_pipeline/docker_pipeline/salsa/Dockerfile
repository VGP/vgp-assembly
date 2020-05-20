FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV SALSA_COMMIT=974589f3302b773dcf0f20c3332fe9daf009fb93
ENV SALSA_DIR=SALSA-2.2-974589

### python
# 2.7
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/2.7.13/Python-2.7.13.tgz && \
    tar xvf Python-2.7.13.tgz && \
    rm /opt/python/Python-2.7.13.tgz && \
    cd Python-2.7.13 && \
    mkdir -p /root/bin/python_2.7.13 && \
    ./configure --with-ensurepip=install --prefix=/root/bin/python_2.7.13 && \
    make install && \
    echo "#%Module" >/root/modules/python/2.7 && \
    echo "prepend-path PATH /root/bin/python_2.7.13/bin" >>/root/modules/python/2.7

### Salsa
WORKDIR /root/tools/salsa
RUN git clone https://github.com/marbl/SALSA.git $SALSA_DIR && \
    cd $SALSA_DIR && \
    git fetch && \
    git checkout $SALSA_COMMIT && \
    make

RUN bash -i -c "module load python/2.7 && \
    which python && \
    which pip && \
    pip install networkx==1.1 numpy "

WORKDIR /root/scripts/salsa/
COPY tmp/* ./
WORKDIR /data
