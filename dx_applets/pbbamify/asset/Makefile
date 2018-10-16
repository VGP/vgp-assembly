all:
	# First change out boost
	sudo apt-get purge -y --auto-remove libboost-dev
	wget http://downloads.sourceforge.net/project/boost/boost/1.67.0/boost_1_67_0.tar.bz2
	tar jxvf boost_1_67_0.tar.bz2
	cd boost_1_67_0/ && ./bootstrap.sh --prefix=/usr/local --with-libraries=all --libdir=/usr/local/lib --includedir=/usr/local/include && sudo ./b2 install
	sudo /bin/bash -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/boost.conf'
	sudo ldconfig

	# Next, install the latest version of cmake
	sudo apt remove -y cmake
	sudo apt purge -y --auto-remove cmake
	mkdir ~/temp && cd ~/temp/ && wget https://cmake.org/files/v3.11/cmake-3.11.4.tar.gz && tar zxvf cmake-3.11.4.tar.gz
	cd ~/temp/cmake-3.11.4/ && ./bootstrap && make -j `nproc` && sudo make install

	# Now download HTSlib
	sudo apt-get install -y liblzma-dev
	wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
	tar jxvf htslib-1.8.tar.bz2
	cd htslib-1.8 && ./configure --prefix=/usr/local --libdir=/usr/local/lib --includedir=/usr/local/include && make -j `nproc` && sudo make install
	sudo /bin/bash -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/htslib.conf'
	sudo ldconfig

	# Now download Google GTest
	git clone https://github.com/google/googletest.git
	mkdir -p googletest/build/ && cd googletest/build/ && cmake ../ && make -j `nproc`

	# Get samtools
	sudo apt-get install -y samtools

	# And grab PacBio pbbam
	git clone https://github.com/PacificBiosciences/pbbam.git
	mkdir -p pbbam/build/ && cd pbbam/build/ && PKG_CONFIG_PATH=/usr/lib/htslib-1.8/lib/pkgconfig/ cmake ../ -DGTEST_SRC_DIR="/home/dnanexus/googletest/" && make -j `nproc`
	sudo cp pbbam/build/bin/* /usr/local/bin/
