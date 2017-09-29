#!/bin/bash

# Install package dependecies
# sudo yum -y install $(cat packages.txt)

# Install virtualenv for deploy a new enviroment
# virtualenv --python=/usr/bin/python2.7  /home/user/Work/Projects/pipeline/.venv
# source /home/user/Work/Projects/pipeline/.venv/bin/activate

# pip install --upgrade pip
# pip install -r modules.txt

# Install sextractor
# echo "saving sextractor package"
# wget -O /home/user/Work/Projects/pipeline/tmp/sextractor.rpm https://www.astromatic.net/download/sextractor/sextractor-2.19.5-1.x86_64.rpm
# sudo rpm -ivh /home/user/Work/Projects/pipeline/tmp/sextractor.rpm
# echo "removing temporal file sextractor.rpm"
# rm /home/user/Work/Projects/pipeline/tmp/sextractor.rpm

# Install scamp
# echo "saving scamp package"
# wget -O /home/user/Work/Projects/pipeline/tmp/scamp.rpm https://www.astromatic.net/download/scamp/scamp-2.0.4-1.x86_64.rpm
# sudo rpm -ivh /home/user/Work/Projects/pipeline/tmp/scamp.rpm
# echo "removing temporal file scamp.rpm"
# rm /home/user/Work/Projects/pipeline/tmp/scamp.rpm

# Install scamp from scratch
# Compile ATLAS/Lapack library
cd /home/user/Work/Projects/pipeline/
mkdir tmp
cd tmp

wget -O atlas3.10.3.tar.bz2 https://downloads.sourceforge.net/project/math-atlas/Stable/3.10.3/atlas3.10.3.tar.bz2?r=&ts=1506698217&use_mirror=10gbps-io
wget -O lapack-3.7.1.tgz http://www.netlib.org/lapack/lapack-3.7.1.tgz

tar -xvf atlas3.10.3.tar.bz2
rm atlas3.10.3.tar.bz2

cd ATLAS
mkdir DONE
cd DONE

../configure -Fa alg -fPIC --with-netlib-lapack-tarfile=../../lapack-3.7.1.tgz
make


# wget -O cdsclient.tar.gz http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz
# tar -xvzf cdsclient.tar.gz
# cd cdsclient-*

# ./configure -prefix=/home/user/Work/Projects/pipeline/.local
# make
# make install

# ./configure --with-atlas-incdir=/home/user/Work/Projects/pipeline/.local/ATLAS/include --with-atlas-libdir=/home/user/Work/Projects/pipeline/.local/ATLAS/lib --without-PLPlot --without-aclient-cgi