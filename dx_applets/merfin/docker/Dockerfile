FROM rocker/r-ver:3.6.1

LABEL com.dnanexus.tool="merfin"



RUN apt update \
  && apt -y install \
    gcc make libc-dev \ 
    zlib1g-dev libbz2-dev liblzma-dev libncurses-dev bzip2 \
  && apt-get install nano \
  && apt-get -y upgrade \
	&& apt-get install -y python3-pip python3-dev build-essential wget bzip2 libz-dev \
	&& ln -s /usr/bin/python3 /usr/bin/python



RUN apt-get install -y git && git clone https://github.com/arangrhie/merfin.git && cd merfin/src && make -j 12



ENTRYPOINT ["/bin/bash", "-c"]
