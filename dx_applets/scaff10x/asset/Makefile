SHELL=/bin/bash -e

all:
	git clone https://github.com/wtsi-hpag/Scaff10X
	cd Scaff10X; bash ./install.sh
	sudo mv Scaff10X/src/scaff-bin /usr/bin/
	sudo mv Scaff10X/src/break10x /usr/bin/
	sudo mv Scaff10X/src/scaff10x /usr/bin/
	sudo mv Scaff10X/src/scaff_reads /usr/bin/