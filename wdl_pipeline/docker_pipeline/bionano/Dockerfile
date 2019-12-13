#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

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

 ### freebayes
 # 1.3.1
 WORKDIR /opt/freebayes
 RUN git clone --recursive git://github.com/ekg/freebayes.git && \
 	cd freebayes && \
 	make && \
 	mkdir /root/bin/freebayes_1.3.1  && \
 	cp bin/freebayes /root/bin/freebayes_1.3.1 && \
 	mkdir /root/modules/freebayes && \
     echo "#%Module" >>/root/modules/freebayes/1.3.1 && \
     echo "append-path PATH /root/bin/freebayes_1.3.1" >>/root/modules/freebayes/1.3.1 && \
     echo "#%Module" >>/root/modules/freebayes/.modulerc && \
     echo "module-version /1.3.1 default" >>/root/modules/freebayes/.modulerc 

# 1.2.0
 WORKDIR /opt/freebayes
 RUN cd freebayes && \
 	git checkout v1.2.0 && \
 	make && \
 	mkdir /root/bin/freebayes_1.2.0  && \
 	cp bin/freebayes /root/bin/freebayes_1.2.0 && \
     echo "#%Module" >>/root/modules/freebayes/1.2.0 && \
     echo "prepend-path PATH /root/bin/freebayes_1.2.0" >>/root/modules/freebayes/1.2.0

### picard 
# 2.9.2
 WORKDIR /opt/picard
 RUN git clone https://github.com/broadinstitute/picard.git && \
     cd picard/ && \
     git checkout 2.9.2 && \
     ./gradlew shadowJar && \
     mv build/libs/picard-2.9.2-SNAPSHOT-all.jar build/libs/picard.jar && \
     mkdir /root/bin/picard_2.9.2 && \
     cp build/libs/picard.jar /root/bin/picard_2.9.2 && \
     mkdir /root/modules/picard && \
     echo "#%Module" >>/root/modules/picard/2.9.2 && \
     echo "setenv PICARD_PATH /root/bin/picard_2.9.2" >> /root/modules/picard/2.9.2

# 2.20.6
 WORKDIR /opt/picard
 RUN cd picard/ && \
     git checkout 2.20.6 && \
     ./gradlew shadowJar && \
     mv build/libs/picard-2.20.6-SNAPSHOT-all.jar build/libs/picard.jar && \
     mkdir /root/bin/picard_2.20.6 && \
     cp build/libs/picard.jar /root/bin/picard_2.20.6 && \
     echo "#%Module" >>/root/modules/picard/2.20.6 && \
     echo "setenv PICARD_PATH /root/bin/picard_2.20.6" >> /root/modules/picard/2.20.6


WORKDIR /root/scripts/bionano/
COPY _docker_bionano.sh tmp/* ./
COPY trf409.linux64 tmp/
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/bionano/_docker_bionano.sh"]