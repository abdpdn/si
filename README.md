## Prerequesites

```
cd /path/to/home
mkdir src
cd src
```

1. Git clone this repository somewhere  : e.g. path-to-home/src/XXZ-Delta-si

2. Install, configure and compile PETSc

2a. Installation

git clone -b release https://gitlab.com/petsc/petsc.git petsc

2b. Configuration: You will need cmake and a recent compiler.
Here I use intel compiler, but gcc should be fine as well. 
You should use the specific tarball of strumpack that you have downloaded. 
Pick a name for the PETSc ARCH, here I use real as PETSc is by default configured to work with real numbers.
Make sure that the option --with-mpi-dir is given the correct location (`which impi`, or alernatively check 
the warning about the environment variable MPI_DIR)
Note : See the bash files for examples.
Note 2: as of February 2024, it was needed to use the development (and not the release) version of petsc and slepc, because of how they fit with Strumpack. 
You can try with the Release of March 2024 (-b release option in git clone)
```
2c. Compilation. Use the command line provided at the end of configuration. 
In my case it was:

```
make PETSC_DIR=/home/jcolbois/src/petsc/ PETSC_ARCH=real all
```
3. Install, configure and compile SLEPc
git clone -b release https://gitlab.com/slepc/slepc
NOTE : see the bash files for examples


4. Boost: find the source of boost somewhere (this is not really essential) and untar
```
cd ..
wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz
tar -xf https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz
```

5. Adapting the cmake configuration options. (Specific to lptsv7)
```
cd path-to-home/src/XXZ-Delta-si
```
5a. Remove any `-xHost` or `-march=native` flags in `CMakeLists.txt` (allows to adapt to the architecture where
the code is executed rather than to that where it is compiled).
5b. If `lptsv7.cmake` is not in `./conf`, copy the one from `olympe`. 
Add a line for `ED_LAPACK_LIB64`:
```
set(ED_LAPACK_LIB64 pthread mkl_lapack95_ilp64 mkl_blas95_ilp64 mkl_core iomp5 mkl_intel_thread mkl_intel_ilp64)
```
Also, make sure that `CMAKE_CXX_COMPILER` and `CMAKE_CC_COMPILER` are set to the correct values.

6. Compiling: adapt docmake.sh as required:
```
module purge
module load intel
module load impi
module load cmake
export CC=mpiicc
export CXX=mpiicpc
export PETSC_DIR=path-to-home/src/petsc/
export PETSC_ARCH=real
export SLEPC_DIR=/path-to-home/src/slepc/
export ED_PETSC_ARCH=$PETSC_ARCH
export ED_SLEPC_ARCH=$PETSC_ARCH
export BOOST_DIR=path-to-home/src/boost_1_80_0/
cmake -DMACHINE=lptsv7 path-to-home/src/XXZ-Delta-si/
```
then 
```
bash docmake.sh
make si
```

(If the make command fail, re-load the modules).

7. LPTSV7 specific subtleties
At the moment when loading the module intel, one obtains
```
$ldd si
...
libquadmath.so.0 => /usr/lib64/libquadmath.so.0 (0x00007f7fce50c000)
...
```
on the login node, and
```
...
libquadmath.so.0 => not found
...
```
on the compute nodes.
As a temporary fix one can copy this library to the build folder and add the local path to the `LD_LIBRARY_PATH`
```
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
``` 
8. Running the code
8a. '`si` admits options which can either be given in command line (e.g. 'spin_si -L 12', ) or in the file `slepc.options`
. Here is an exemple:
```
-L 16
-Sz 0
-disorder 5.0
-seed 3
-Delta 0.08

-pbc 1
-eps_nev 100

-measure_correlations
-measure_entanglement
-measure_local

-st_type sinvert

-st_pc_factor_mat_solver_package strumpack
-st_pc_type lu
-mat_strumpack_verbose
-mat_strumpack_colperm 0
```
Command lines options take precedence over the ones in this file.

8b. `srun` options on lptsv7:
```
srun -N nodes --ntasks=ntasks --ntasks-per-node=n -c ncpu-per-task ./si
```
