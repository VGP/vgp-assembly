all:
	# download HTSlib
	sudo apt-get install -y liblzma-dev
	wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
	tar jxvf htslib-1.9.tar.bz2
	cd htslib-1.9 && ./configure --prefix=/usr/local --libdir=/usr/local/lib --includedir=/usr/local/include && make -j `nproc` && sudo make install
	sudo /bin/bash -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/htslib.conf'
	sudo ldconfig
	
	# download BCFtools
	wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
	tar jxvf bcftools-1.9.tar.bz2
	cd bcftools-1.9 && ./configure --prefix=/usr/local --libdir=/usr/local/lib --includedir=/usr/local/include && make -j `nproc` && sudo make install
	
	# download samtools
	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	tar jxvf samtools-1.9.tar.bz2
	cd samtools-1.9 && ./configure --prefix=/usr/local --libdir=/usr/local/lib --includedir=/usr/local/include && make -j `nproc` && sudo make install
