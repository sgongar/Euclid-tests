#!/bin/bash
# An small script for libraries intallation
# 
# Google's Shell Style Guide


function upgrade_system {
  sudo yum -y install $(cat packages.txt)
}


function install_virtualenv {
  Install virtualenv for deploy a new enviroment
  virtualenv --python=/usr/bin/python2.7  /home/user/Work/Projects/pipeline/.venv
  source /home/user/Work/Projects/pipeline/.venv/bin/activate
}


function update_pip {
  pip install --upgrade pip
  pip install -r modules.txt
}


function install_atlas {
  ATLAS_URL="https://downloads.sourceforge.net/project/math-atlas/Stable/3.10.3/atlas3.10.3.tar.bz2?r=&ts=1506698217&use_mirror=10gbps-io"
  LAPACK_URL="http://www.netlib.org/lapack/lapack-3.7.1.tgz"

  wget -O atlas.tar.bz2 $ATLAS_URL
  wget -O lapack.tgz $LAPACK_URL

  tar -xf atlas.tar.bz2
  rm atlas.tar.bz2

  cd ATLAS
  if [ ! -d DONE/ ]; then
    mkdir DONE
  fi
  cd DONE/

  mkdir /home/user/Work/Projects/pipeline/.local/ATLAS
  ../configure -Fa alg -fPIC --with-netlib-lapack-tarfile=../../lapack.tgz --prefix=/home/user/Work/Projects/pipeline/.local/ATLAS
  make

  # Get into lib folder and make shared libraries
  cd lib
  make shared
  cd ../

  make install
}


function install_cdsclient {
  cdsclient_url="http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz"
  wget -O cdsclient.tar.gz $cdsclient_url

  if [ ! -d cdsclient/ ]; then
    mkdir cdsclient/
  fi

  tar -xzf cdsclient.tar.gz -C cdsclient --strip-components=1
  cd cdsclient

  # Create a dir for cdsclient installation
  cdsclient_dir="/home/user/Work/Projects/pipeline/.local/cdsclient"

  if [ ! -d $cdsclient_dir ]; then
    mkdir $cdsclient_dir
  fi

  # Configure
  ./configure -prefix=/home/user/Work/Projects/pipeline/.local/cdsclient

  # Compile them
  make

  # Perfom an installation to local folder
  make install
}


function install_sextractor {
  SEXTRACTOR_URL="https://www.astromatic.net/download/sextractor/sextractor-2.19.5.tar.gz"
  wget -O sextractor.tar.gz $SEXTRACTOR_URL

  if [ ! -d sextractor/ ]; then
    mkdir sextractor/
  fi

  tar -xzf sextractor.tar.gz -C sextractor/ --strip-components=1
  cd sextractor

  ATLAS_INCLUDE_DIR="/home/user/Work/Projects/pipeline/.local/ATLAS/include"
  ATLAS_LIB_DIR="/home/user/Work/Projects/pipeline/.local/ATLAS/lib"
  SEXTRACTOR_BIN_DIR="/home/user/Work/Projects/pipeline/.local"
  ./configure --with-atlas-incdir=$ATLAS_INCLUDE_DIR \
  --with-atlas-libdir=$ATLAS_LIB_DIR --prefix=$SEXTRACTOR_BIN_DIR

  make
  make install
}


function install_scamp {
  SCAMP_URL="https://www.astromatic.net/download/scamp/scamp-2.0.4.tar.gz"
  wget -O scamp.tar.gz $SCAMP_URL

  if [ ! -d scamp/ ]; then
    mkdir scamp/
  fi

  tar -xzf scamp.tar.gz -C scamp/ --strip-components=1
  cd scamp

  ATLAS_INCLUDE_DIR="/home/user/Work/Projects/pipeline/.local/ATLAS/include"
  ATLAS_LIB_DIR="/home/user/Work/Projects/pipeline/.local/ATLAS/lib"
  CDSCLIENT_BIN_DIR="/home/user/Work/Projects/pipeline/.local/cdsclient/bin"
  SCAMP_BIN_DIR="/home/user/Work/Projects/pipeline/.local"

  ./configure --with-atlas-incdir=$ATLAS_INCLUDE_DIR \
  --with-atlas-libdir=$ATLAS_LIB_DIR --with-cdsclient-dir=$CDSCLIENT_BIN_DIR \
  --prefix=$SCAMP_BIN_DIR

  make
  make install
}


function copy_files {
  cp -r /media/sf_Euclid-tests/pipeline/* /home/user/Work/Projects/pipeline/* 
}


function main {
  TMP_DIR="/home/user/Work/Projects/pipeline/tmp"
  LOCAL_DIR="/home/user/Work/Projects/pipeline/.local/"

  # Checking directories
  if [ ! -d "$TMP_DIR" ]; then
    mkdir $TMP_DIR
  fi

  if [ ! -d "$LOCAL_DIR" ]; then
    mkdir $LOCAL_DIR
  fi

  # Install scamp from scratch
  # Compile ATLAS/Lapack library
  cd $TMP_DIR 

  upgrade_system

  install_virtualenv

  update_pip

  install_atlas

  cd ../../
  rm -rf *

  install_cdsclient
  cd ../

  install_sextractor
  cd ../

  install_scamp
  cd ../
  rm -rf scamp*
}


if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  main "$@"
fi
