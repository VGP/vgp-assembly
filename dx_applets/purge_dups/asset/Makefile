all:
	curl -L https://github.com/lh3/minimap2/releases/download/v2.16/minimap2-2.16_x64-linux.tar.bz2 && mv minimap2-2.16_x64-linux.tar.bz2 / && cd / && tar --no-same-owner -jxvf minimap2-2.16_x64-linux.tar.bz2 
	git clone https://github.com/dfguan/purge_dups.git && mv purge_dups / && cd /purge_dups/src && make
	git clone https://github.com/dfguan/runner.git && mv purge_dups / && cd /runner && python3 setup.py install --user
