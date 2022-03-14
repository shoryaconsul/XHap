git clone https://github.com/vibansal/HapCUT2.git
cd HapCUT2
git clone https://github.com/samtools/htslib.git
cd htslib
git submodule update --init --recursive
autoreconf -i
./configure --disable-bz2 --disable-lzma
make
cd ..
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$(pwd)/htslib"
make
