FROM centos:7.5.1804
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV SOLVE_CODE_FILENAME=Solve3.3_10252018.tar.gz

RUN yum update -y && \
    yum install -y perl perl-Data-Dumper epel-release R python3 && \
    echo -e "#!/bin/bash\necho \"All prereqs are installed in PATH. Skipping 'module' call: '\$0 \$@'\"\n" >>/usr/bin/module && \
    chmod +x /usr/bin/module

WORKDIR /opt/bionano/
COPY tmp/$SOLVE_CODE_FILENAME ./
RUN tar xvf $SOLVE_CODE_FILENAME && \
    rm $SOLVE_CODE_FILENAME && \
    cd Solve3.3_10252018

WORKDIR /root/scripts/bionano/
COPY tmp/*.sh ./

WORKDIR /root/scripts/bionano/trimNs/
COPY tmp/trimNs/* ./

WORKDIR /root/config/bionano
COPY tmp/*.xml ./

WORKDIR /data
