static char help[] =
    "Shift-Invert ED for XXZ chain with Fibonnaci fields \n (C) Fabien Alet "
    "2017. \n\n"
    "The command line options are:\n"
    "  -L <L>, where <L> = size of the chain [default 6]\n"
    "  ( .... Description not complete ....) \n"
    "\n";

// #define PETSC_DESIRE_COMPLEX
// #define PETSC_USE_COMPLEX
// #include <complex>
#include <omp.h>
#include <slepceps.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <random>
#ifdef USE_MKL
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#else
typedef struct __MKL_Complex16 {
  double real;
  double imag;
} MKL_Complex16;
typedef int MKL_INT;
typedef MKL_Complex16 lapack_complex_double;

#include <cblas.h>
#include <lapacke.h>
#endif

double PI = acos(-1.);
using namespace std;

#include "Spin/Spin_basis.h"

#include "Spin/Spin_parameters.h"

#include "Spin/SpinOneHalfXXZ_disorder.h"

#include "Spin/Spin_observable.h"

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv) {
  cout.precision(20);
  SlepcInitialize(&argc, &argv, "slepc.options", help);

  /**** Init parallel work ****/
  // For parallelization on node (openMP)
  int ENV_NUM_THREADS = mkl_get_max_threads();  /// Get value of OMP_NUM_THREADS
  omp_set_num_threads(
      ENV_NUM_THREADS);  // Deactivate OpenMP for PETSc and SLEPc
#ifdef USE_MKL
  mkl_set_num_threads(
      ENV_NUM_THREADS);  // Deactivate OpenMP for PETSc and SLEPc
#endif

  // For parallelization between nodes (MPI)
  int myrank, mpisize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  /***** Petsc / Slepc data structure ***/
  Mat H;
  EPS eps;
  EPSType type;
  PetscErrorCode ierr;
  PetscReal error, tol, re, im;

  int number_of_states = 2;
  double Sz = 0;
  ierr = PetscOptionsGetReal(NULL, NULL, "-Sz", &Sz, NULL);
  PetscCall(ierr);

  Parameters myparameters(myrank);
  int L = myparameters.L;
  basis mybasis(L, Sz, myrank, number_of_states);
  //  mybasis.check_basis(); // ??
  int nconf = mybasis.total_number_of_confs;

  if (myrank == 0) {
    std::cout << "# L = " << L << " number of states=" << nconf << std::endl;
  }
  observable myobservable(&mybasis, ENV_NUM_THREADS);

  PetscInt Istart, Iend;
  /************************** Hamiltonian **********************************/
  {
    std::vector<int> local_block_sizes(mpisize, nconf / mpisize);
    for (size_t i = 0; i < nconf % mpisize; i++)
      local_block_sizes[i]++;  // distribute evenly
    Istart = 0;
    for (size_t i = 0; i < myrank; i++) Istart += local_block_sizes[i];
    Iend = Istart + local_block_sizes[myrank];

    Hamiltonian myHamiltonian(&mybasis, &H);
    myHamiltonian.get_parameters();

    myparameters.string_from_basis = mybasis.string_names;
    myparameters.string_from_H = myHamiltonian.string_names;

    //	 myobservable.get_local_field_x(myHamiltonian.local_field_x);
    // myobservable.get_local_field_z(myHamiltonian.local_field_z);
    /*myHamiltonian.init_matrix_sizes(Istart, Iend);
    ierr = MatSetUp(H);
    PetscCall(ierr);
    myHamiltonian.create_matrix(Istart, Iend);

    if (myrank == 0){ 
      myHamiltonian.print_local_fields();
    }*/

  }

  /****  Matrix assembly ****/
  /*MatSetOption(H, MAT_SYMMETRIC, PETSC_TRUE);
  MatSetOption(H, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);
  if (myrank == 0) std::cout << "# Assembly... " << std::flush;
  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  PetscCall(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  PetscCall(ierr);
  if (myrank == 0) {
    std::cout << "Hamiltonian matrix assembly done." << std::endl;
  }

  // MatView(H,PETSC_VIEWER_STDOUT_WORLD);
  PetscBool output_matrix=PETSC_FALSE;
  PetscOptionsGetBool(NULL, NULL,"-output_matrix",&output_matrix,NULL);

    if (output_matrix) {
      PetscViewer    view;
      PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Matrix.mtx",&view);
      PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATRIXMARKET);
      MatView(H,view);
      PetscViewerPopFormat(view);
      PetscViewerDestroy(&view);
      exit(0);
    }
  */


  // TO UNCOMMENT
  if (myrank == 0) {
    std::cout << "#L= " << L << " ,";
    std::cout << " number of states= " << nconf;
    std::cout << "\n";
  } 
  
  PetscScalar Eminc, Emaxc;
  Eminc = 0; Emaxc = 1; 
  
  /*bool do_ev = 1;
  if (do_ev) {
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
    PetscCall(ierr);
    ierr = EPSSetOperators(eps, H, NULL);
    PetscCall(ierr);
    ierr = EPSSetProblemType(eps, EPS_HEP);
    PetscCall(ierr);
    EPSSetDimensions(eps, 1, PETSC_DECIDE, PETSC_DECIDE);

    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    ierr = EPSSolve(eps);
    PetscCall(ierr);
    EPSGetEigenvalue(eps, 0, &Eminc, NULL);
    double Emin = PetscRealPart(Eminc);
    if (0 == myrank) std::cout << "Emin of H = " << Emin << "\n";
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
    ierr = EPSSolve(eps);
    PetscCall(ierr);
    EPSGetEigenvalue(eps, 0, &Emaxc, NULL);
    double Emax = PetscRealPart(Emaxc);
    if (0 == myrank) std::cout << "Emax of H = " << Emax << "\n";
  }*/

  /************* Measurement on eigenvectors *********/
  // Initialize file names

  // TO UNCOMMENT
  // Vec xr, xi;
  // MatCreateVecs(H, PETSC_NULLPTR, &xr);
  // MatCreateVecs(H, PETSC_NULLPTR, &xi);
  double renorm_target = 0.5;
  /*
  std::vector<double> targets;
  if (myparameters.target_infinite_temperature) {
    PetscScalar E_infinite_temperature;
    ierr = MatGetTrace(H, &E_infinite_temperature);
    PetscCall(ierr);
    E_infinite_temperature /= mybasis.total_number_of_confs;
    if (0 == myrank)
      std::cout << "E(beta = 0) = " << E_infinite_temperature << " E_renorm = "
                << (E_infinite_temperature - Eminc) / (Emaxc - Eminc)
                << std::endl;
    targets.push_back((E_infinite_temperature - Eminc) / (Emaxc - Eminc));
  } else {
    std::cout << "Given a list of targets... " << std::endl;
    targets = myparameters.targets;
  }*/

  // new setup
  /*
  double target = renorm_target * (Emaxc - Eminc) + Eminc;
  EPS eps2;
  ST st;
  PetscInt nconv = 0;
   
  int trys = 0;
  bool diagfailed = true;
  */

  // TO UNCOMMENT

  int nconv = 1;
  std::cout << "# Going through converged eigenstates.... nconv =  " << nconv << std::endl;
    

  ofstream entout;
  ofstream locout;
  ofstream partout;
  ofstream corrout;
  ofstream tcorrout;


  myparameters.init_filenames_eigenstate(entout, locout, partout,corrout, tcorrout,
                                             renorm_target);


  std::vector<double> energies;
  // std::vector<double> rgap;

  PetscScalar Er, Ei;
  Vec Vec_local;
  std::vector<double> Siz(L, 0.);
  std::vector< std::vector<double> > Gij(L);
  for (int k=0;k<L;++k) { Gij[k].resize(L,0.);}
  std::vector< std::vector<double> > Tij(L);
  for (int k=0;k<L;++k) { Tij[k].resize(L,0.);}


  PetscScalar *state;
  state = (PetscScalar *)calloc(nconf, sizeof(PetscScalar));
  PetscInt seed = 15101976;
  PetscOptionsGetInt(NULL, NULL, "-seed", &seed, NULL);
  PetscInt num_states = 100;
  PetscOptionsGetInt(NULL, NULL, "-num_states", &num_states, NULL);
  PetscCall(ierr);

    // random generator

  boost::mt19937 random_generator(seed);
  boost::random::uniform_real_distribution<double> box_dist(0, 1);
  
    // create the states
  if (myrank == 0) {
    std::cout << "# Creating " << num_states << " random states." << std::endl;
  }

  for (int r = 0; r < num_states; r++) {
    // create a random Sz0 
    double aa = 0;
    double bb = 0;
    double val;

      
      // TO UNCOMMENT
    double norm = 0;
      

    for (int i = 0; i < nconf; i++) {
      while (aa == 0) {
        aa = box_dist(random_generator);
      }
      bb = box_dist(random_generator);
      val = sqrt(-2.0 * log(aa)) * cos(2 * acos(-1.) * bb);
      while (aa == 0) {
        aa = box_dist(random_generator);
      }
      bb = box_dist(random_generator);
      state[i] = (PetscScalar)val;
      norm += val * val;
    }

    // // normalize
    for (int i = 0; i < nconf; i++) {
      state[i] /= sqrt(norm);
    }


    //ierr = EPSGetEigenpair(eps2, i, &Er, &Ei, xr, xi);
    //PetscCall(ierr);
       
    energies.push_back(0.5);
    if (myparameters.write_wf) {
      std::stringstream filename;
      filename << "Eigenvector" << r << "."
             << myparameters.string_from_basis
             << myparameters.string_from_H << ".dat";
      std::string str(filename.str());
          // write vector 
          /*
          PetscViewer viewer;
          PetscViewerBinaryOpen(PETSC_COMM_WORLD, str.c_str(), FILE_MODE_WRITE,
                                &viewer);
          VecView(xr, viewer);
          */
    }
        
        /*if (mpisize == 1) {
          VecCreateSeq(PETSC_COMM_SELF, nconf, &Vec_local);
          ierr = VecCopy(xr, Vec_local);
        }  // only 1 mpi proc
        else {
          VecScatter ctx;
          VecScatterCreateToZero(xr, &ctx, &Vec_local);
          VecScatterBegin(ctx, xr, Vec_local, INSERT_VALUES, SCATTER_FORWARD);
          VecScatterEnd(ctx, xr, Vec_local, INSERT_VALUES, SCATTER_FORWARD);
          VecScatterDestroy(&ctx);
        }*/
        
        
    if (myrank == 0) {
          // PetscScalar *state;
          // VecGetArray(Vec_local, &state);


      if (myparameters.measure_all_part_entropy) {
        double qmin = myparameters.qmin;
        std::vector<double> qs;
        double mymin = qmin;
        double current_q;
        if (qmin < 1) {
          current_q = qmin;
          for (int p = 0; p < 9; ++p) {
            qs.push_back(current_q);
            current_q += (1. - qmin) / 9.;
          }
          mymin = 1;
        }
        current_q = mymin;
        double dq = (myparameters.qmax - mymin) / (myparameters.Nq - 1);
        for (int nq = 0; nq < myparameters.Nq - 1; ++nq) {
            qs.push_back(current_q);
            current_q += dq;
        }
        qs.push_back(current_q);
        double real_Nq = qs.size();
        std::vector<double> entropies(real_Nq, 0.);
        entropies = myobservable.all_part_entropy(
        state, myparameters.Nq, myparameters.qmin, myparameters.qmax);

            for (int nq = 0; nq < real_Nq; ++nq) {
              partout << "S[q=" << qs[nq] << "] " << entropies[nq] << "\n";
            }
            partout << "S[q=infty]= " << myobservable.Smax(state) << "\n";
      } else {
            double PE = myobservable.part_entropy(state, 1);
            // parout << PetscRealPart(Er) << " " << PE << "\n";}
            partout << PE << "\n";
      }

          if (myparameters.measure_local) {
            myobservable.compute_local_magnetization(state);
            Siz = myobservable.sz_local;
            for (int r = 0; r < L; ++r) {
              locout << r << " " << Siz[r] << endl;
            }
          }

          if (myparameters.measure_correlations) {
            Gij =myobservable.get_two_points_connected_correlation(state);
            for (int r = 0; r < L; ++r) {
              for (int s = r; s < L; ++s) {
	              corrout << r << " " << " " << s << " " << Gij[r][s] << endl;
              } 
            }

          }

          if (myparameters.measure_transverse_correlations){
            // get two points transverse correlations
            Tij = myobservable.get_SpSm_correlation(state);
            for (int r = 0; r < L; ++r) {
              for (int s = r; s < L; ++s) {
                tcorrout << r << " " << " " << s << " " << Tij[r][s] << endl;
              } 
            }
          }

          /*
          if (myparameters.measure_entanglement) {
            myobservable.compute_entanglement_spectrum(state);
            double S1 = myobservable.entang_entropy(1);
            entout << S1 << "\n";
          }
          */
          // VecRestoreArray(Vec_local, &state);
        
    }

           // VecDestroy(&Vec_local);
  }  

  if (myrank == 0) {
    std::cout << "# Done with eigenstates." << std::endl;
    entout.close();
    locout.close();
    partout.close();
    corrout.close();
    tcorrout.close();
  }

  if (myrank == 0) {

    cout << "# Energies " << endl;
    ofstream enout;
    ofstream rgapout;
    //double target = renorm_target * (Emaxc - Eminc) + Eminc;
    double target = 0.5;
    double renorm_target = 0.5;
    myparameters.init_filenames_energy(enout, rgapout, renorm_target);
    // as a header add info about the min/max energies
    double Eminc = 0; double Emaxc = 1;
    enout << "# (Emin, Emax) = " << Eminc << "," << Emaxc << endl;
    enout << "# Etarget = " << target << endl;
    enout << "# Erenorm_target = " << renorm_target << endl;
    
    for (int i = 0; i < nconv*num_states; i++) {
      cout << energies[i] << "\n";
      enout << energies[i] << endl;
    }
    
    std::sort(energies.begin(), energies.end());

    enout.close();
    rgapout.close();
  }
  
  //  }
  SlepcFinalize();
  return 0;
}
