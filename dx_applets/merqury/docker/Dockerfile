FROM rocker/r-ver:3.6.1

LABEL com.dnanexus.tool="merqury"


RUN apt-get update -qq && apt-get -y --no-install-recommends install \
	wget \
	libssl-dev \  
	libcurl4-openssl-dev \
	&& rm -rf /var/lib/apt/lists/* \
  	&& R -e "install.packages('argparse', repos ='http://cran.us.r-project.org', dependencies=TRUE)" \
  	&& R -e "install.packages('ggplot2', repos ='http://cran.us.r-project.org', dependencies=TRUE)" \
  	&& R -e "install.packages('scales', repos ='http://cran.us.r-project.org', dependencies=TRUE)"


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
  && make install \
  && apt-get install nano \
  && apt-get -y upgrade \
	&& apt-get install -y python3-pip python3-dev build-essential wget bzip2 libz-dev \
	&& ln -s /usr/bin/python3 /usr/bin/python

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz -O /tmp/bedtools.tar.gz && \
	tar zxvf /tmp/bedtools.tar.gz -C /opt/ && rm /tmp/bedtools.tar.gz && \
	cd /opt/bedtools2 && \
	make

ENV PATH="/opt/bedtools2/bin/:${PATH}"


# Setup ENV variables
ENV \
  IGVTOOLS_HOME=/opt/IGVTools \
  IGVTOOLS_VERSION=2.3.98

# Install IGVTools
RUN \
  wget --quiet -O igvtools_${IGVTOOLS_VERSION}.zip \
    http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_${IGVTOOLS_VERSION}.zip \
  && unzip igvtools_${IGVTOOLS_VERSION}.zip \
  && rm igvtools_${IGVTOOLS_VERSION}.zip \
  && mv IGVTools $IGVTOOLS_HOME

ENV PATH="${IGVTOOLS_HOME}:${PATH}"

# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw


RUN wget https://github.com/marbl/meryl/releases/download/v1.0/meryl-1.0.Linux-amd64.tar.xz -O /tmp/meryl-1.0.Linux-amd64.tar.xz && \
	tar -xJf /tmp/meryl-1.0.Linux-amd64.tar.xz -C /opt/ && rm /tmp/meryl-1.0.Linux-amd64.tar.xz 
ENV PATH="/opt/meryl-1.0/Linux-amd64/bin:${PATH}"


# Install OpenJDK-8 & Fix certificate issues
RUN apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;


# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN apt-get install -y git && git clone https://github.com/marbl/merqury.git && cd merqury && git checkout tags/v1.1
ENV PATH="$PWD:${PATH}"
ENV MERQURY="/merqury"


ENTRYPOINT ["/bin/bash", "-c"]
