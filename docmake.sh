module purge
module load intel
module load impi
module load cmake
export CC=mpiicc
export CXX=mpiicpc
export PETSC_DIR=/home/jcolbois/src/petsc/
export PETSC_ARCH=real2
export SLEPC_DIR=/home/jcolbois/src/slepc/
export ED_PETSC_ARCH=$PETSC_ARCH
export ED_SLEPC_ARCH=$PETSC_ARCH
export BOOST_DIR=/home/jcolbois/src/boost_1_80_0/
cmake -DMACHINE=lptsv7 /home/jcolbois/src/XXZ-Delta-si/
