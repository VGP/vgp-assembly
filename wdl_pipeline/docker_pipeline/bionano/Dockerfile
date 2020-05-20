FROM centos:7.5.1804
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV SOLVE_CODE_FILENAME=Solve3.3_10252018.tar.gz

RUN yum update -y && \
    yum install -y perl perl-Data-Dumper epel-release R python3 && \
    pip3 install biopython && \
    echo -e "#!/bin/bash\necho \"Prereqs are installed in docker image PATH. Skipping 'module' call: '\$0 \$@'\"\n" >>/usr/bin/module && \
    chmod +x /usr/bin/module

WORKDIR /opt/bionano/
COPY tmp/$SOLVE_CODE_FILENAME ./
RUN tar xvf $SOLVE_CODE_FILENAME && \
    rm $SOLVE_CODE_FILENAME && \
    cd Solve3.3_10252018
RUN echo "Bionano Genomics has agreed to provide the licensed Bionano Solve software to enable " >/opt/bionano/README && \
    echo "the VGP consortium to package the VGP pipeline in a container. The Bionano Solve " >>/opt/bionano/README && \
    echo "software in the VGP pipeline may not be the latest version, may not have the most " >>/opt/bionano/README && \
    echo "recent security patches and may be configured in a way unsupported by Bionano. The " >>/opt/bionano/README && \
    echo "VGP fully assumes the maintenance and support of this VGP pipeline. Please reach out " >>/opt/bionano/README && \
    echo "to tpesout@ucsc.edu with issues. For the latest supported Bionano software, " >>/opt/bionano/README && \
    echo "visit Bionano Support (https://bionanogenomics.com/support/software-downloads/)." >>/opt/bionano/README

WORKDIR /root/scripts/bionano/
COPY tmp/*.sh ./

WORKDIR /root/scripts/bionano/trimNs/
COPY tmp/trimNs/* ./

WORKDIR /root/config/bionano
COPY tmp/*.xml ./

WORKDIR /data
