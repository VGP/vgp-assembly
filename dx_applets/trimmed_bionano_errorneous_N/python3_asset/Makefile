SHELL=/bin/bash -e -x -o pipefail
all:
	
	sudo apt-get update
		
	sudo apt-get install -y libssl-dev openssl
	sudo wget https://www.python.org/ftp/python/3.4.0/Python-3.4.0.tgz
	sudo tar xzvf Python-3.4.0.tgz
	pushd Python-3.4.0 && sudo ./configure && sudo make && sudo make install && popd;
	python3.4 -m pip install biopython

