
cd petsc_git
module purge
module load intel
module load impi
module load cmake
export CC=mpiicc
export CXX=mpiicpc

make PETSC_DIR=/home/jcolbois/src/petsc_git PETSC_ARCH=real all

make PETSC_DIR=/home/jcolbois/src/petsc_git PETSC_ARCH=real check
