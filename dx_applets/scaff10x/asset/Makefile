SHELL=/bin/bash -e

all:
	git clone https://github.com/wtsi-hpag/Scaff10X
	cd Scaff10X; git checkout tags/v4.1
	cd Scaff10X; bash ./install.sh;
	# for some reason scaff_FilePreProcess won't compile succesfully
	cd Scaff10X/src/scaff-bin; wget --no-check-certificate https://github.com/wtsi-hpag/Scaff10X/releases/download/v4.1/scaff_FilePreProcess -O scaff_FilePreProcess;
	sudo chmod 777 Scaff10X/src/scaff-bin/scaff_FilePreProcess
	sudo mv Scaff10X/src/scaff-bin /usr/bin/
	sudo mv Scaff10X/src/break10x /usr/bin/
	sudo mv Scaff10X/src/scaff10x /usr/bin/
	sudo mv Scaff10X/src/scaff_reads /usr/bin/