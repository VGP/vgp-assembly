FROM ubuntu:18.04

RUN apt-get update && \
    apt-get install -y wget  zlib1g-dev && \
    apt install build-essential -y --no-install-recommends


RUN wget https://github.com/chhylp123/hifiasm/archive/v0.11.tar.gz && tar -zxvf v0.11.tar.gz && cd /hifiasm-0.11/ && make

 