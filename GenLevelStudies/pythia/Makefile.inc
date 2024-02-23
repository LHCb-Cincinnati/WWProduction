# PYTHIA configuration file.
# Generated on Fri Feb 23 09:59:47 EST 2024 with the user supplied options:
# --with-root-bin=/usr/local/Cellar/root/6.30.04/bin
# --with-python-config=/usr/local/Cellar/python@3.11/3.11.7_2/bin/python3.11-config

# Install directory prefixes.
PREFIX_BIN=/Applications/pythia8310/bin
PREFIX_INCLUDE=/Applications/pythia8310/include
PREFIX_LIB=/Applications/pythia8310/lib
PREFIX_SHARE=/Applications/pythia8310/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
LIB_SUFFIX=.dylib
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/usr/local/Cellar/root/6.30.04/bin/
ROOT_INCLUDE=-I/usr/local/Cellar/root/6.30.04/include/root
ROOT_LIB=-L/usr/local/Cellar/root/6.30.04/lib/root -Wl,-rpath,/usr/local/Cellar/root/6.30.04/lib/root -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=true
PYTHON_CONFIG=/usr/local/Cellar/python@3.11/3.11.7_2/bin/python3.11-config
PYTHON_BIN=
PYTHON_INCLUDE=-I/usr/local/Cellar/python@3.11/3.11.7_2/Frameworks/Python.framework/Versions/3.11/include/python3.11 -I/usr/local/Cellar/python@3.11/3.11.7_2/Frameworks/Python.framework/Versions/3.11/include/python3.11
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
