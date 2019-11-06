#TODO FROM tpesout/modules:4.3.0--c999abdc958fe45c735dca093555b5413aecb678
FROM tpesout/modules:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

############################################################
# Finished 
# # need to create a versioned folder within gnuplot even though we will always default to it
# WORKDIR /opt/gnuplot
# RUN wget https://sourceforge.net/projects/gnuplot/files/gnuplot/5.2.7/gnuplot-5.2.7.tar.gz && \
# 	tar xvf gnuplot-5.2.7.tar.gz && \
# 	cd gnuplot-5.2.7 && \
# 	./configure && \
# 	make && \
# 	make install && \
# 	mkdir /root/bin/gnuplot_5.2.7  && \
# 	cp -r . /root/bin/gnuplot_5.2.7 && \
# 	rm -r /opt/gnuplot && \
# 	mkdir /root/modules/gnuplot && \
# 	echo "#%Module" >>/root/modules/gnuplot/5.2.7 && \
#     echo "append-path PATH /root/bin/gnuplot_5.2.7" >>/root/modules/gnuplot/5.2.7 && \
#     echo "#%Module" >>/root/modules/gnuplot/.modulerc && \
#     echo "module-version /5.2.7 default" >>/root/modules/gnuplot/.modulerc

# WORKDIR /opt/bwa
# RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 && \
# 	tar xvf bwa-0.7.17.tar.bz2 && \
# 	cd bwa-0.7.17 && \
# 	# ./configure && \
# 	make && \
# 	# make install && \
# 	mkdir /root/bin/bwa_0.7.17  && \
# 	cp -r . /root/bin/bwa_0.7.17 && \
# 	rm -r /opt/bwa && \
# 	mkdir /root/modules/bwa && \
# 	echo "#%Module" >>/root/modules/bwa/0.7.17 && \
#     echo "append-path PATH /root/bin/bwa_0.7.17" >>/root/modules/bwa/0.7.17 && \
#     echo "#%Module" >>/root/modules/bwa/.modulerc && \
#     echo "module-version /0.7.17 default" >>/root/modules/bwa/.modulerc 

# WORKDIR /opt/mash
# RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar && \
# 	tar xvf mash-Linux64-v2.2.tar && \
# 	cd mash-Linux64-v2.2 && \
# 	# ./bootstrap.sh && \
# 	# ./configure && \
# 	# make && \
# 	# make install && \
# 	mkdir /root/bin/mash_Linux64_v2.2  && \
# 	cp -r . /root/bin/mash_Linux64_v2.2 && \
# 	rm -r /opt/mash && \
# 	mkdir /root/modules/mash && \
# 	echo "#%Module" >>/root/modules/mash/2.2 && \
#     echo "append-path PATH /root/bin/mash_Linux64_v2.2" >>/root/modules/mash/2.2 && \
#     echo "#%Module" >>/root/modules/mash/.modulerc && \
#     echo "module-version /2.2 default" >>/root/modules/mash/.modulerc 

# WORKDIR /opt/mummer
# RUN wget https://github.com/mummer4/mummer/releases/download/v3.9.4alpha/mummer-3.9.4alpha.tar.gz && \
# 	tar xvf mummer-3.9.4alpha.tar.gz && \
# 	cd mummer-3.9.4alpha && \
# 	./configure && \
# 	make && \
# 	make install && \
# 	ldconfig && \
# 	mkdir /root/bin/mummer_3.9.4alpha  && \
# 	cp -r . /root/bin/mummer_3.9.4alpha && \
# 	rm -r /opt/mummer && \
# 	mkdir /root/modules/mummer && \
# 	echo "#%Module" >>/root/modules/mummer/3.9.4 && \
#     echo "append-path PATH /root/bin/mummer_3.9.4alpha" >>/root/modules/mummer/3.9.4 && \
#     echo "#%Module" >>/root/modules/mummer/.modulerc && \
#     echo "module-version /3.9.4 default" >>/root/modules/mummer/.modulerc 

# WORKDIR /opt/ngmlr
# RUN wget https://github.com/philres/ngmlr/releases/download/v0.2.6/ngmlr-0.2.6-beta-linux-x86_64.tar.gz && \
# 	tar xvf ngmlr-0.2.6-beta-linux-x86_64.tar.gz && \
# 	cd ngmlr-0.2.6 && \
# 	# ./configure && \
# 	# make && \
# 	# make install && \
# 	mkdir /root/bin/ngmlr_0.2.6_beta_linux_x86_64  && \
# 	cp -r . /root/bin/ngmlr_0.2.6_beta_linux_x86_64 && \
# 	rm -r /opt/ngmlr && \
# 	mkdir /root/modules/ngmlr && \
# 	echo "#%Module" >>/root/modules/ngmlr/0.2.6 && \
#     echo "append-path PATH /root/bin/ngmlr_0.2.6_beta_linux_x86_64" >>/root/modules/ngmlr/0.2.6 && \
#     echo "#%Module" >>/root/modules/ngmlr/.modulerc && \
#     echo "module-version /0.2.6 default" >>/root/modules/ngmlr/.modulerc

# WORKDIR /opt/minimap2
# RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 && \
# 	tar xvf minimap2-2.17.tar.bz2 && \
# 	cd minimap2-2.17 && \
# 	make && \
# 	# make install && \
# 	mkdir /root/bin/minimap2_2.17  && \
# 	cp -r minimap2 /root/bin/minimap2_2.17 && \
# 	rm -r /opt/minimap2 && \
# 	# mkdir /root/modules/minimap2 && \
# 	echo "#%Module" >>/root/modules/minimap2/2.17 && \
#     echo "append-path PATH /root/bin/minimap2_2.17" >>/root/modules/minimap2/2.17
#     # echo "#%Module" >>/root/modules/minimap2/.modulerc && \
#     # echo "module-version /5.2.7 default" >>/root/modules/minimap2/.modulerc

############################################################
# working on 

# Multiple version installation must be handled specially
# WORKDIR /opt/python3
# RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
#     tar xvf Python-3.6.0.tgz && \
#     cd Python-3.6.0 && \
#     ./configure && \
#     make && \
#     mkdir /root/bin/python_3.6.0 && \
#     cp python /root/bin/python_3.6.0 && \
#     rm -r /opt/python3 && \
#     mkdir /root/modules/python && \
#     echo "#%Module" >>/root/modules/python/3.6 && \
#     echo "append-path PATH /root/bin/python_3.6.0" >>/root/modules/python/3.6 && \
#     echo "#%Module" >>/root/modules/python/.modulerc && \
#     echo "module-version /3.6 default" >>/root/modules/python/.modulerc

# WORKDIR /opt/python2
# RUN wget https://www.python.org/ftp/python/2.7/Python-2.7.tgz && \
#     tar xvf Python-2.7.tgz && \
#     cd Python-2.7 && \
#     ./configure && \
#     make && \
#     mkdir /root/bin/python_2.7 && \
#     cp python /root/bin/python_3.6.0 && \
#     rm -r /opt/python3 && \
#     mkdir /root/modules/python && \
#     echo "#%Module" >>/root/modules/python/3.6 && \
#     echo "append-path PATH /root/bin/python_3.6.0" >>/root/modules/python/3.6 && \
#     echo "#%Module" >>/root/modules/python/.modulerc && \
#     echo "module-version /3.6 default" >>/root/modules/python/.modulerc

# WORKDIR /opt/perl
# RUN wget https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/perl/5.18.2-2ubuntu1/perl_5.18.2.orig.tar.bz2 && \
# 	tar xvf perl_5.18.2.orig.tar.bz2 && \
# 	cd perl-5.18.2 && \
# 	./Configure -des -Dprefix=/root/modules/perl && \
# 	make && \
# 	make install && \
# 	mkdir /root/bin/perl_5.18.2  && \
# 	cp -r . /root/bin/perl_5.18.2 && \
# 	rm -r /opt/perl && \
# 	mkdir /root/modules/perl && \
# 	echo "#%Module" >>/root/modules/perl/5.18.2 && \
#     echo "append-path PATH /root/bin/perl_5.18.2" >>/root/modules/perl/5.18.2 && \
#     echo "#%Module" >>/root/modules/perl/.modulerc && \
#     echo "module-version /5.18.2 default" >>/root/modules/perl/.modulerc


############################################################
# template
# WORKDIR /opt/python3
# RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
#     tar xvf Python-3.6.0.tgz && \
#     cd Python-3.6.0 && \
#     ./configure && \
#     make && \
#     mkdir /root/bin/python_3.6.0 && \
#     cp python /root/bin/python_3.6.0 && \
#     rm -r /opt/python3 && \
#     mkdir /root/modules/python && \
#     echo "#%Module" >>/root/modules/python/3.6 && \
#     echo "append-path PATH /root/bin/python_3.6.0" >>/root/modules/python/3.6 && \
#     echo "#%Module" >>/root/modules/python/.modulerc && \
#     echo "module-version /3.6 default" >>/root/modules/python/.modulerc

# TODO later

# WORKDIR /opt/repeatmasker
# RUN wget http://www.repeatmasker.org/RepeatMasker-open-4-0-9-p2.tar.gz && \
# 	tar xvf RepeatMasker-open-4-0-9-p2.tar.gz && \
# 	cd RepeatMasker && \
# 	# ./configure && \
# 	#make && \
# 	#make install && \
# 	mkdir /root/bin/RepeatMasker_open_4_0_9_p2  && \
# 	cp -r . /root/bin/RepeatMasker_open_4_0_9_p2 && \
# 	rm -r /opt/repeatmasker && \
# 	mkdir /root/modules/repeatmasker && \
# 	echo "#%Module" >>/root/modules/repeatmasker/4.0.9 && \
#     echo "append-path PATH /root/bin/RepeatMasker_open_4_0_9_p2" >>/root/modules/repeatmasker/4.0.9 && \
#     echo "#%Module" >>/root/modules/repeatmasker/.modulerc && \
#     echo "module-version /4.0.9 default" >>/root/modules/repeatmasker/.modulerc 

# WORKDIR /opt/picard
# RUN wget https://github.com/broadinstitute/picard/archive/2.9.2.tar.gz && \
#     tar xvf 2.9.2.tar.gz && \
#     cd picard-2.9.2 && \
#	  ./gradlew shadowJar
#     ./configure && \
#     make && \
#     mkdir /root/bin/python_3.6.0 && \
#     cp python /root/bin/python_3.6.0 && \
#     rm -r /opt/python3 && \
#     mkdir /root/modules/python && \
#     echo "#%Module" >>/root/modules/python/3.6 && \
#     echo "append-path PATH /root/bin/python_3.6.0" >>/root/modules/python/3.6 && \
#     echo "#%Module" >>/root/modules/python/.modulerc && \
#     echo "module-version /3.6 default" >>/root/modules/python/.modulerc

WORKDIR /root/scripts/mashmap/
COPY _docker_mashmap.sh tmp/*.sh ./
WORKDIR /data
ENTRYPOINT ["bash", "-i", "/root/scripts/mashmap/_docker_mashmap.sh"]
