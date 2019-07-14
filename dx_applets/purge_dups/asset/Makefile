all:
	wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && mv minimap2-2.17_x64-linux.tar.bz2 / && cd / && tar --no-same-owner -jxvf minimap2-2.17_x64-linux.tar.bz2 
	git clone https://github.com/dfguan/purge_dups.git && mv purge_dups / && cd /purge_dups/src && make
	apt-get install python3-pip -y
	pip install setuptools==3.0 
	git clone https://github.com/dfguan/runner.git && mv runner / && cd /runner && python3 setup.py install --user
