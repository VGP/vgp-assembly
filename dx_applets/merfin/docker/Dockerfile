FROM rocker/r-ver:3.6.1

LABEL com.dnanexus.tool="merfin"


RUN apt update \
  && apt -y install \
    gcc make libc-dev \ 
    zlib1g-dev libbz2-dev liblzma-dev libncurses-dev bzip2 \
  && apt-get install nano \
  && apt-get -y upgrade \
  && apt-get install -y tabix \
	&& apt-get install -y python3-pip python3-dev build-essential wget bzip2 libz-dev \
	&& ln -s /usr/bin/python3 /usr/bin/python

RUN apt update && apt upgrade && apt install curl 

RUN apt-get install -y git && git clone https://github.com/arangrhie/merfin.git && cd merfin/src && make -j 12

RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar -vxjf htslib-1.9.tar.bz2 && cd htslib-1.9 && make

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar -vxjf bcftools-1.9.tar.bz2 && cd bcftools-1.9 && make

ENV PATH="/bcftools-1.9/:${PATH}"

ENTRYPOINT ["/bin/bash", "-c"]
