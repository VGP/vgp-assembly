#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# need to create a versioned folder within gnuplot even though we will always default to it
WORKDIR /opt/gnuplot
RUN wget https://sourceforge.net/projects/gnuplot/files/gnuplot/5.2.7/gnuplot-5.2.7.tar.gz && \
	tar xvf gnuplot-5.2.7.tar.gz && \
	cd gnuplot-5.2.7 && \
	./configure && \
	make && \
	make install && \
	mkdir /root/bin/gnuplot_5.2.7  && \
	cp -r . /root/bin/gnuplot_5.2.7 && \
	rm -r /opt/gnuplot && \
	mkdir /root/modules/gnuplot && \
	echo "#%Module" >>/root/modules/gnuplot/5.2.7 && \
    echo "append-path PATH /root/bin/gnuplot_5.2.7" >>/root/modules/gnuplot/5.2.7 && \
    echo "#%Module" >>/root/modules/gnuplot/.modulerc && \
    echo "module-version /5.2.7 default" >>/root/modules/gnuplot/.modulerc

WORKDIR /opt/bwa
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 && \
	tar xvf bwa-0.7.17.tar.bz2 && \
	cd bwa-0.7.17 && \
	make && \
	mkdir /root/bin/bwa_0.7.17  && \
	cp -r . /root/bin/bwa_0.7.17 && \
	rm -r /opt/bwa && \
	mkdir /root/modules/bwa && \
	echo "#%Module" >>/root/modules/bwa/0.7.17 && \
    echo "append-path PATH /root/bin/bwa_0.7.17" >>/root/modules/bwa/0.7.17 && \
    echo "#%Module" >>/root/modules/bwa/.modulerc && \
    echo "module-version /0.7.17 default" >>/root/modules/bwa/.modulerc 

WORKDIR /opt/mash
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar && \
	tar xvf mash-Linux64-v2.2.tar && \
	cd mash-Linux64-v2.2 && \
	mkdir /root/bin/mash_Linux64_v2.2  && \
	cp -r . /root/bin/mash_Linux64_v2.2 && \
	rm -r /opt/mash && \
	mkdir /root/modules/mash && \
	echo "#%Module" >>/root/modules/mash/2.2 && \
    echo "append-path PATH /root/bin/mash_Linux64_v2.2" >>/root/modules/mash/2.2 && \
    echo "#%Module" >>/root/modules/mash/.modulerc && \
    echo "module-version /2.2 default" >>/root/modules/mash/.modulerc 

WORKDIR /opt/mummer
RUN wget https://github.com/mummer4/mummer/releases/download/v3.9.4alpha/mummer-3.9.4alpha.tar.gz && \
	tar xvf mummer-3.9.4alpha.tar.gz && \
	cd mummer-3.9.4alpha && \
	./configure && \
	make && \
	make install && \
	ldconfig && \
	mkdir /root/bin/mummer_3.9.4alpha  && \
	cp -r . /root/bin/mummer_3.9.4alpha && \
	rm -r /opt/mummer && \
	mkdir /root/modules/mummer && \
	echo "#%Module" >>/root/modules/mummer/3.9.4 && \
    echo "append-path PATH /root/bin/mummer_3.9.4alpha" >>/root/modules/mummer/3.9.4 && \
    echo "#%Module" >>/root/modules/mummer/.modulerc && \
    echo "module-version /3.9.4 default" >>/root/modules/mummer/.modulerc 

WORKDIR /opt/ngmlr
RUN wget https://github.com/philres/ngmlr/releases/download/v0.2.6/ngmlr-0.2.6-beta-linux-x86_64.tar.gz && \
	tar xvf ngmlr-0.2.6-beta-linux-x86_64.tar.gz && \
	cd ngmlr-0.2.6 && \
	mkdir /root/bin/ngmlr_0.2.6_beta_linux_x86_64  && \
	cp -r . /root/bin/ngmlr_0.2.6_beta_linux_x86_64 && \
	rm -r /opt/ngmlr && \
	mkdir /root/modules/ngmlr && \
	echo "#%Module" >>/root/modules/ngmlr/0.2.6 && \
    echo "append-path PATH /root/bin/ngmlr_0.2.6_beta_linux_x86_64" >>/root/modules/ngmlr/0.2.6 && \
    echo "#%Module" >>/root/modules/ngmlr/.modulerc && \
    echo "module-version /0.2.6 default" >>/root/modules/ngmlr/.modulerc

WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 && \
	tar xvf minimap2-2.17.tar.bz2 && \
	cd minimap2-2.17 && \
	make && \
	mkdir /root/bin/minimap2_2.17  && \
	cp -r minimap2 /root/bin/minimap2_2.17 && \
	rm -r /opt/minimap2 && \
	echo "#%Module" >>/root/modules/minimap2/2.17 && \
    echo "append-path PATH /root/bin/minimap2_2.17" >>/root/modules/minimap2/2.17

WORKDIR /root/scripts/mashmap/
COPY _docker_mashmap.sh tmp/*.sh ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/mashmap/_docker_mashmap.sh"]
