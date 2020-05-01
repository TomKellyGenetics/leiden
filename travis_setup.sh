#!/bin/bash

# install python
if [[ $PYENV == "py27"]]; then
    version = "2.7"
    if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
    elif [[ $TRAVIS_OS_NAME == "osx" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh
   fi
fi
if [[ $PYENV == "py35"]]; then
    version "3.5"
    if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh-O miniconda.sh
    elif [[ $TRAVIS_OS_NAME == "osx" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh-O miniconda.sh
   fi
fi
if [[ $PYENV == "py36"]]; then
    version "3.6"
    if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    elif [[ $TRAVIS_OS_NAME == "osx" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
   fi
fi
if [[ $PYENV == "py37"]]; then
    version "3.7"
    if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -O miniconda.sh
    elif [[ $TRAVIS_OS_NAME == "osx" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-MacOSX-x86_64.sh -O miniconda.sh
   fi
fi
if [[ $PYENV == "py38"]]; then
    version "3.8"
    if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-py38_4.8.2-Linux-x86_64.sh -O miniconda.sh
    elif [[ $TRAVIS_OS_NAME == "osx" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-py38_4.8.2-MacOSX-x86_64.sh -O miniconda.sh
   fi
fi

bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
export RETICULATE_PYTHON="$HOME/miniconda/bin/python"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update --user -q conda
conda info -a

conda init $0
conda create -y -n r-reticulate r-base igraph python=$version
conda activate r-reticulate
pip install --user --upgrade pip
pip install --user igraph leidenalg
conda install -c conda-forge leidenalg

