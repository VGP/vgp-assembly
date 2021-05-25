FROM ubuntu:18.04

RUN apt-get update && \
    apt-get install -y wget  zlib1g-dev && \
    apt install build-essential -y --no-install-recommends


RUN wget https://github.com/chhylp123/hifiasm/archive/0.15.2.tar.gz && tar -zxvf 0.15.2.tar.gz && cd /hifiasm-0.15.2/ && make

 