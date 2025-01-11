
cd slepc_git
module purge
module load intel
module load impi
module load cmake

export PETSC_DIR=/home/jcolbois/src/petsc_git
export PETSC_ARCH=real
export SLEPC_DIR=/home/jcolbois/src/slepc_git


make SLEPC_DIR=/home/jcolbois/src/slepc_git PETSC_DIR=/home/jcolbois/src/petsc_git PETSC_ARCH=real
make SLEPC_DIR=/home/jcolbois/src/slepc_git PETSC_DIR=/home/jcolbois/src/petsc_git check
