
cd petsc_git
module purge
module load intel
module load impi

source /opt/intel/oneapi/mpi/2021.5.1/env/vars.sh

module load cmake


./configure --force --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" --with-blas-lapack-dir=$MKLROOT --with-debugging=0 --with-errorchecking=0 --download-scalapack --download-parmetis --download-metis --download-strumpack --with-openmp --with-precision=double --with-make-np=1 MPI-DIR=$MPI_DIR PETSC_ARCH=real-main


