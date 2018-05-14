#!/bin/bash
# An small script for libraries installation
#
# *GNU Terry Pratchett*

# Should be already installed
# sudo yum install cmake
# sudo yum install blas-devel lapack-devel could be useful?

# Elements 5.2.1
# Download it!
make configure
make install
# move to directory

# Alexandria 2.9
https://gitlab.euclid-sgs.uk/hdegaude/Alexandria/tree/8f6d8331704b314b0e221860f2e05a9f564a8fa6
make configure
make
make install

# Lapack
# Copy make.inc.example to make.inc

# Sextractor
# Requires:
# - levmar
# - opencv / from official repository
# - opencv-devel / from official repository
# - yaml-cpp / from official repository
# - yaml-cpp-devel / from official repository

# sudo yum install glibc


make configure
