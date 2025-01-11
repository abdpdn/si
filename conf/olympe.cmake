# For mumps double precision
set(ED_USE_MKL 1)
set(ED_MKLROOT       $ENV{MKLROOT})
set(ED_PETSC_DIR     $ENV{PETSC_DIR})
set(ED_SLEPC_DIR     $ENV{SLEPC_DIR})
set(ED_PETSC_ARCH     $ENV{PETSC_ARCH})
set(ED_SLEPC_ARCH     $ENV{SLEPC_ARCH})


set(ED_BOOST_ROOT    $ENV{BOOST_DIR})

set(ED_MPI_HOME      $ENV{I_MPI_ROOT})

set(ED_LAPACK_LIB_DIR $ENV{MKLROOT}/lib/intel64/)
#set(ED_LAPACK_LIB     mkl_scalapack_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread)
set(ED_LAPACK_LIB     mkl_lapack95_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread)
set(ED_LAPACK_LIB64 pthread mkl_lapack95_ilp64 mkl_blas95_ilp64 mkl_core iomp5 mkl_intel_thread mkl_intel_ilp64)

#set(ED_LAPACK_LIB     mkl_scalapack_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core mkl_blacs_openmpi_lp64 iomp5 pthread)

set(CMAKE_CXX_COMPILER "mpiicpc")
set(CMAKE_CC_COMPILER "mpiicc")
