SHELL=/bin/bash -e -x -o pipefail

all:
	# Install minimap2
	curl -L https://github.com/lh3/minimap2/releases/download/v2.9/minimap2-2.9_x64-linux.tar.bz2 | tar jxvf -; ( sudo mv ./minimap2*/minimap2 /usr/bin/ )
	# Install biobambam2 
	sudo apt-get install libtool -y
	wget https://github.com/gt1/libmaus2/archive/2.0.455-release-20180319165218.tar.gz -O -|tar zx; ( cd libmaus2-2.0.455-release-20180319165218; libtoolize --force; aclocal; autoreconf -i -f; ./configure; make; make install; rm -rf ../libmaus2-2.0.455-release-20180319165218 )
	# libmaus2 is installed under /usr/local/
	wget https://github.com/gt1/biobambam2/archive/2.0.87-release-20180301132713.tar.gz -O -|tar zx; ( cd biobambam2-2.0.87-release-20180301132713; autoreconf -i -f; ./configure --with-libmaus2="/usr/local/lib" --prefix="/usr/local"; make; make install; rm -rf ../biobambam2-2.0.87-release-20180301132713 )
	# Install sambamba
	wget https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2 -O - | tar xj; ( sudo mv sambamba /usr/local/bin/ )
