SHELL=/bin/bash -e -o pipefail -x

all:
	wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && mv minimap2-2.17_x64-linux.tar.bz2 / && cd / && tar --no-same-owner -jxvf minimap2-2.17_x64-linux.tar.bz2 
	#pip install --upgrade pip 
	pip install setuptools==3.0 && apt-get install python3-pip -y  && git clone https://github.com/dfguan/runner.git && mv runner /therunner && cd /therunner && python3 setup.py install 
	git clone https://github.com/dfguan/purge_dups.git && cd purge_dups && git checkout debug && cd .. && mv purge_dups / && cd /purge_dups/src && make
