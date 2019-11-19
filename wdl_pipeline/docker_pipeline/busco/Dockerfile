FROM tpesout/vgp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

### BUSCO
# 3.0.2
WORKDIR /root/tools/BUSCO
RUN wget https://gitlab.com/ezlab/busco/-/archive/3.0.2/busco-3.0.2.tar.gz && \
    tar xvf busco-3.0.2.tar.gz && \
    mv busco-3.0.2 busco/ && \
    rm busco-3.0.2.tar.gz && \
    cd busco && \
    bash -i -c "module load python/3.6 && \
    python setup.py install"
COPY config.ini /root/tools/BUSCO/busco/config/config.ini
# vertebrata dataset
WORKDIR /root/tools/BUSCO/dataset/
RUN wget http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz && \
    tar xvf vertebrata_odb9.tar.gz && \
    rm vertebrata_odb9.tar.gz

### Augustus
# 3.3.1
WORKDIR /root/tools/Augustus/
RUN wget https://github.com/Gaius-Augustus/Augustus/archive/v3.3.1-tag1.tar.gz && \
    tar xvf v3.3.1-tag1.tar.gz && \
    mv Augustus-3.3.1-tag1/ Augustus-3.3.1/ && \
    rm v3.3.1-tag1.tar.gz && \
    cd Augustus-3.3.1/ && \
    mkdir tooldir && \
    ln -s /opt/samtools/samtools-1.9/ tooldir/samtools  && \
    ln -s /opt/bcftools/bcftools-1.9/ tooldir/bcftools && \
    ln -s /opt/htslib/htslib-1.9/ tooldir/htslib && \
    TOOLDIR=/root/tools/Augustus/Augustus-3.3.1/tooldir make

WORKDIR /root/scripts/busco/
COPY tmp/*.sh ./
WORKDIR /data
