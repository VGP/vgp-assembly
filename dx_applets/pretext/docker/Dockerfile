FROM continuumio/miniconda

LABEL com.dnanexus.tool="pretext-suite"

RUN conda install -c bioconda pretext-suite

ADD http://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 /tmp/apps/samtools.tar.bz2

RUN apt update \
  && apt -y install \
    gcc make libc-dev \ 
    zlib1g-dev libbz2-dev liblzma-dev libncurses-dev bzip2 \
  && cd /tmp/apps \ 
  && mkdir samtools \
  && tar xjvf samtools.tar.bz2 -C samtools --strip-components=1 \
  && cd samtools \
  && ./configure \
  && make \
  && make install
