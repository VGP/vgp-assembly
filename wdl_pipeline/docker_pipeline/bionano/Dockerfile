#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


# ### python
# 3.5.0
WORKDIR /opt/python3
RUN wget https://www.python.org/ftp/python/3.5.0/Python-3.5.0.tgz && \
    tar xvf Python-3.5.0.tgz && \
    rm -r /opt/python3/Python-3.5.0.tgz && \
    cd Python-3.5.0 && \
    ./configure && \
    make && \
    mkdir /root/bin/python_3.5.0 && \
    cp python /root/bin/python_3.5.0 && \
    # mkdir /root/modules/python && \
    echo "#%Module" >>/root/modules/python/3.5 && \
    echo "prepend-path PATH /root/bin/python_3.5.0" >>/root/modules/python/3.5

### python
# 2.7
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/2.7/Python-2.7.tgz && \
    tar xvf Python-2.7.tgz && \
    rm -r /opt/python/Python-2.7.tgz && \
    cd Python-2.7 && \
    ./configure && \
    make && \
    mkdir /root/bin/python_2.7 && \
    cp python /root/bin/python_2.7 && \
    echo "#%Module" >>/root/modules/python/2.7 && \
    echo "prepend-path PATH /root/bin/python_2.7" >>/root/modules/python/2.7




### perl
# 5.18.2
WORKDIR /opt/perl
RUN wget https://www.cpan.org/src/5.0/perl-5.18.2.tar.gz && \
	tar -xvf  perl-5.18.2.tar.gz && \
	rm -r /opt/perl/perl-5.18.2.tar.gz && \
	cd perl-5.18.2 && \
	./Configure -des -Dprefix=/opt/perl/perl-5.18.2 && \
	make && \
	make install && \
	mkdir /root/bin/perl_5.18.2 && \
	cp perl /root/bin/perl_5.18.2 && \
	mkdir /root/modules/perl && \
	echo "#%Module" >>/root/modules/perl/5.18.2 && \
    echo "prepend-path PATH /root/bin/perl_5.18.2" >>/root/modules/perl/5.18.2

WORKDIR /root/scripts/bionano/
COPY _docker_bionano.sh tmp/* ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/bionano/_docker_bionano.sh"]
