#!/bin/bash

# Install package dependecies
sudo yum -y install $(cat packages.txt)

# Install virtualenv for deploy a new enviroment
virtualenv --python=/usr/bin/python2.7  /home/user/Work/Projects/pipeline/.venv
source /home/user/Work/Projects/pipeline/.venv/bin/activate

pip install --upgrade pip
pip install -r modules.txt

# Install sextractor
echo "saving sextractor package"
wget -O /home/user/Work/Projects/pipeline/tmp/sextractor.rpm https://www.astromatic.net/download/sextractor/sextractor-2.19.5-1.x86_64.rpm
sudo rpm -ivh /home/user/Work/Projects/pipeline/tmp/sextractor.rpm
echo "removing temporal file sextractor.rpm"
rm /home/user/Work/Projects/pipeline/tmp/sextractor.rpm

# Create local directory to deploy libraries and binaries
mkdir /home/user/Work/Projects/pipeline/.local/

# Install scamp from scratch
# Compile ATLAS/Lapack library
cd /home/user/Work/Projects/pipeline/tmp

ATLAS_URL="https://downloads.sourceforge.net/project/math-atlas/Stable/3.10.3/atlas3.10.3.tar.bz2?r=&ts=1506698217&use_mirror=10gbps-io"
LAPACK_URL="http://www.netlib.org/lapack/lapack-3.7.1.tgz"

wget -O atlas.tar.bz2 $ATLAS_URL
wget -O lapack.tgz $LAPACK_URL

tar -xf atlas.tar.bz2
rm atlas.tar.bz2

cd ATLAS
mkdir DONE
cd DONE

../configure -Fa alg -fPIC --with-netlib-lapack-tarfile=../../lapack.tgz
make

# Get into lib folder and make shared libraries
cd lib
make shared
cd ../

awk '{ if (NR == 3) print "DESTDIR=/home/user/Work/Projects/pipeline/.local/ATLAS"; else print $0}' Makefile > output_file.txt
mv output_file.txt Makefile

make install

cd ../../
rm -rf *

wget -O cdsclient.tar.gz http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz
tar -xzf cdsclient.tar.gz
cd cdsclient-*

mkdir /home/user/Work/Projects/pipeline/.local/cdsclient

./configure -prefix=/home/user/Work/Projects/pipeline/.local/cdsclient
make
make install

cd ../
rm -rf cdsclient*

wget -O scamp.tar.gz https://www.astromatic.net/download/scamp/scamp-2.0.4.tar.gz
tar -xzf scamp.tar.gz

cd scamp*

ATLAS_INCLUDE_DIR="/home/user/Work/Projects/pipeline/.local/ATLAS/include"
ATLAS_LIB_DIR="/home/user/Work/Projects/pipeline/.local/ATLAS/lib"
CDSCLIENT_BIN_DIR="/home/user/Work/Projects/pipeline/.local/cdsclient/bin"
SCAMP_BIN_DIR="/home/user/Work/Projects/pipeline/.local"
./configure --with-atlas-incdir=$ATLAS_INCLUDE_DIR --with-atlas-libdir=$ATLAS_LIB_DIR --with-cdsclient-dir=$CDSCLIENT_BIN_DIR --prefix=$SCAMP_BIN_DIR

make
make install

cd ../
rm -rf scamp*
