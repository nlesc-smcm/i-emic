/******************************x****************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

/***************************************************
 * Put together/extended by                        *
 * T.E. Mulder, Utrecht University 2014/15/16/17   *
 * contact: t.e.mulder@uu.nl                       *
 ***************************************************/

// A namespace with utils from hymls and trios that I need

#ifndef UTILS_H
#define UTILS_H

#include "GlobalDefinitions.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <fstream>
#include <vector>
#include <memory>

class Combined_MultiVec;

template<typename Vector>

class ComplexVector;

class Epetra_Comm;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_MultiVector;
class Epetra_CrsMatrix;

namespace EpetraExt {
    class HDF5;
}

namespace TRIOS {
    class Domain;
}

// Perhaps some lapack is nice
extern "C" void dscal_(int* N, double *DA, double *X, int *INCX);

extern "C" void daxpy_(int* N,    double *DA, double *X,
                       int *INCX, double *Y,  int *INCY);

extern "C" double ddot_(int *N, double *X, int *INCX, double *Y, int *INCY);


namespace Utils
{
    //! A simple struct to store CRS info
    struct CRSMat
    {
        std::vector<double> co;
        std::vector<int>    jco;
        std::vector<int>    beg;
    };

    //! We need both a distributed and a global version of the land mask, so
    //! we use the following struct.
    struct MaskStruct
    {
        //! distributed
        Teuchos::RCP<Epetra_IntVector> local;

        //! global
        std::shared_ptr<std::vector<int> > global;

        //! global borderless
        std::shared_ptr<std::vector<int> > global_borderless;

        //! global surface (borderless)
        std::shared_ptr<std::vector<int> > global_surface;

        //! mask name
        std::string label;

        //! dimension in x, y, z-direction
        int n, m, l;
    };

    //! Clone an Epetra_Vector
    Teuchos::RCP<Epetra_Vector> clone(
        Teuchos::RCP<const Epetra_Vector> const &vec);

    //! Clone a Combined_MultiVec
    std::shared_ptr<Combined_MultiVec> clone(
        std::shared_ptr<const Combined_MultiVec> const &vec);

    //! Obtain 2-norm of std::vector<double>
    double norm(std::vector<double> &vec);

    //! Obtain first 2-norm of Epetra_MultiVector
    double norm(Epetra_MultiVector &vec);

    //! Obtain first 2-norm of vectorPtr, more convenient interface
    template<typename VectorPtr>
    double norm(VectorPtr const &vec)
    {
        std::vector<double> norm(vec->NumVectors(), 0.0);
        CHECK_ZERO(vec->Norm2(&norm[0]));
        return norm[0];
    }

    //! Obtain first 2-norm of Combined_MultiVec
    double norm(Combined_MultiVec &vec);

    //! Obtain first inf-norm of Epetra_MultiVector
    double normInf(Epetra_MultiVector &vec);

    //! Obtain first inf-norm of Combined_MultiVec
    double normInf(Combined_MultiVec &vec);

    //! Obtain first inf-norm of vectorPtr, more convenient interface
    template<typename VectorPtr>
    double normInf(VectorPtr const &vec)
    {
        std::vector<double> norm(vec->NumVectors());
        CHECK_ZERO(vec->NormInf(&norm[0]));
        return norm[0];
    }

    //! Return the norm for a specific field specified by <unknown>
    //! (1-6 for an ocean vector)
    double normOfField(Teuchos::RCP<Epetra_MultiVector> vec,
                       Teuchos::RCP<TRIOS::Domain> domain, int unknown);

    //! Compute column sums of matrix mat and put them in vector sums
    void colSums(Epetra_CrsMatrix const &mat, Epetra_Vector &sums);

    //! Update std::vector<double>, result is stored in B
    //! B = scalarA*A+scalarB*B
    void update(double scalarA, std::vector<double> &A,
                double scalarB, std::vector<double> &B);

    //! Print vectorPtr
    template<typename VectorPtr>
    void print(VectorPtr const &vec, std::string const &fname)
    {
        // This does not necessarily behave well in parallel,
        // use proc-dependent filenames or a gather on the vector.
        std::ofstream file;
        file.open(fname);
        vec->Print(file);
        file.close();
    }

    //! Print std vector
    void print(std::vector<double> const &vec, std::string const &fname);

    //! Print surface mask to INFO and file fname
    void printSurfaceMask(MaskStruct const &mask);
    
    //! Print surface mask to INFO and file fname
    void printSurfaceMask(std::shared_ptr<std::vector<int> > mask,
                          std::string const &fname,
                          int nDim);

    //! Compute dot product of two vectors
    template<typename VectorPtr>
    double dot(VectorPtr const &vec1,
               VectorPtr const &vec2)
    {
        double result = 0;
        CHECK_ZERO(vec1->Dot(*vec2, &result));
        return result;
    }

    //! Compute dot product of two vectors
    double dot(std::vector<double> &vec1, std::vector<double> &vec2);

    //! Compute sum of a vector
    double sum(std::vector<double> &vec);

    int nnz(Teuchos::RCP<Epetra_Vector> vec, double t = 1e-12);

    //! The parameters in dominantPars will overwrite parameters in originalPars
    //! that have the same name.
    void overwriteParameters(Teuchos::RCP<Teuchos::ParameterList> originalPars,
                             Teuchos::RCP<Teuchos::ParameterList> dominantPars);

    //! Update a parameterlist from a file
    Teuchos::RCP<Teuchos::ParameterList> obtainParams(std::string const &str,
                                                      std::string const &name);

    //! Update a sublist from a file
    void obtainParams(Teuchos::RCP<Teuchos::ParameterList> pars,
                      std::string const &str,
                      std::string const &name);

    //! Hashing an Epetra_MultiVector
    size_t hash(Teuchos::RCP<Epetra_MultiVector> vec);

    //! Save/load
    void save(Teuchos::RCP<const Epetra_MultiVector> vec, std::string const &filename);
    void load(Teuchos::RCP<Epetra_MultiVector> vec, std::string const &filename);

    //! Save/load a combined multivector to several hdf5 binaries
    void save(std::shared_ptr<const Combined_MultiVec> vec, std::string const &filename);
    void load(std::shared_ptr<Combined_MultiVec> vec, std::string const &filename);
    
    void save(Combined_MultiVec const &vec, std::string const &filename);

    //----------------------------------------------------------------------
    void saveEigenvectors(std::vector<ComplexVector<Combined_MultiVec> > const &eigvs,
                          std::vector<std::complex<double> > const &alpha,
                          std::vector<std::complex<double> > const &beta,
                          std::string const &filename);

    //----------------------------------------------------------------------
    void saveEigenvectors(std::vector<ComplexVector<Epetra_Vector> > const &eigvs,
                          std::vector<std::complex<double> > const &alpha,
                          std::vector<std::complex<double> > const &beta,
                          std::string const &filename);

    //----------------------------------------------------------------------
    void saveEigenvalues(EpetraExt::HDF5 &HDF5,
                         std::vector<std::complex<double> > const &alpha,
                         std::vector<std::complex<double> > const &beta,
                         int numEigs);

    //-----------------------------------------------------------------------------
    template <typename EigenSolver>
    void saveEigenvectors(EigenSolver eigSolver, std::string const &filename)
    {
        auto eigvs = eigSolver->getEigenVectors();
        auto alpha = eigSolver->getAlpha();
        auto beta  = eigSolver->getBeta();
        
        saveEigenvectors(eigvs, alpha, beta, filename);
    }

    //!------------------------------------------------------------------
    //! getVector wrapper to select view or copy (for Teuchos::RCP's)
    Teuchos::RCP<Epetra_Vector> getVector(char mode,
                                          Teuchos::RCP<Epetra_Vector> vec);

    //!------------------------------------------------------------------
    //! getVector wrapper to select view or copy (for shared_ptr's)
    std::shared_ptr<std::vector<double> > getVector
    (char mode, std::shared_ptr<std::vector<double> > vec);

    //!------------------------------------------------------------------
    //! Assemble local or global crs struct into an Epetra_CrsMatrix.
    //! This handles both 0 and 1-based CRS structs.
    void assembleCRS(Teuchos::RCP<Epetra_CrsMatrix> mat,
                     CRSMat const &crs, int const maxnnz,
                     Teuchos::RCP<TRIOS::Domain> domain = Teuchos::null);
    
    //------------------------------------------------------------------
    // Original Utils from Jonas
    //------------------------------------------------------------------
    
    //------------------------------------------------------------------
    //! split a cartesian box into nparts subdomains in a somehow 'good' way.
    //! This function just tells you how many subdomains there should be in
    //! every direction. Returns 1 if no splitting was found.
    int SplitBox(int nx, int ny, int nz,
                 int nparts, int& ndx, int& ndy, int& ndz,
                 int sx = 1, int sy = 1, int sz = 1);

    //! converts linear index to cartesian subscripts
    void ind2sub(int nx, int ny, int nz, int dof,
                 int idx, int& i, int& j, int& k, int& var);
    //! converts linear index to cartesian subscripts
    void ind2sub(int nx, int ny, int nz, int idx, int& i, int& j, int& k);
    //! converts cartesian subscripts to linear index
    int sub2ind(int nx, int ny, int nz, int dof, int i, int j, int k, int var);

    //! create a map for a 3D problem with multiple degrees of freedom per node

    /*! creates a distributed non-overlapping map to be used for constructing/solving
      linear problems with one or more dof/node. This is similar to the Galeri map
      Cartesian3D, except that we allow several dof/node.
      For 2D problems simply set nz = 1.

      This function can also be used if the number of grid cells nx*ny*nz
      is less than the number of processes, in which case any proc with PID
      >= nx*ny*nz gets 0 elements in the map.

      It is possible to artifically restrict the number of procs that get any nodes to
      maxNumActiveProcs by setting the flag to something >0.
      The default -1 means all in comm get nodes (if there are sufficient nodes).
    */
    Teuchos::RCP<Epetra_Map> CreateMap(int nx, int ny, int nz,
                                       int dof, int indexbase,
                                       const Epetra_Comm& comm,
                                       int maxNumActiveProcs=-1);

    //! a general map creating function for a 3D grid

    /* creates a map of the form (i0:i1,j0:j1,k0:k1),
       column-major access (i is fastest).
       Small letters indicate the bounds of the
       subdomain, capital letters the bounds of the
       global domain. The indices of subdomains are
       allowed to overlap                          */
    Teuchos::RCP<Epetra_Map> CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
                                       int I0, int I1, int J0, int J1, int K0, int K1,
                                       int dof,
                                       const Epetra_Comm& comm);

    //! creates a map of the form (i0:i1,j0:j1,k0:k1),
    //! column-major access (i is fastest).
    //! Small letters indicate the bounds of the
    //! subdomain, capital letters the bounds of the
    //! global domain. The indices of subdomains are
    //! allowed to overlap                          */
    Teuchos::RCP<Epetra_Map> CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
                                       int I0, int I1, int J0, int J1, int K0, int K1,
                                       const Epetra_Comm& comm);

    //! extract a map with all indices of the original one between min_gid and max_gid
    //! (including the bounds)
    Teuchos::RCP<Epetra_Map> ExtractRange(const Epetra_Map& map, int min_gid, int max_gid);

    //! extract a map with 1 variable from a map with dof variables per node. 'var'
    //! is the variable to be extracted, i.e. 0, 1,...dof-1
    Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, int dof, int var);

    //! extract a map for two variables from a map with dof variables. 'var'
    //! are the variables to be extracted, i.e. {0,1}, {4,5} etc.
    Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, int dof, const int var[2]);

    //! extract a map for nvars variables from a map with dof variables. 'var'
    //! are the variables to be extracted.
    Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, int dof,
                                          const int *var, int nvars);

    //! given a map and an array indicating wether each node of the map is to be
    //! discarded (true) or not (false), this function creates a new map with the
    //! discarded entries removed.
    Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, const bool* discard);

    //! Extract a map based on a list of global indices
    Teuchos::RCP<Epetra_BlockMap> CreateSubMap(const Epetra_BlockMap& map,
                                               std::vector<int> const &indices);

    //! Given an Epetra_MultiVector and a list of global indices we return an Epetra_MultiVector with
    //! the same distribution but restricted to the values at the supplied indices.
    Teuchos::RCP<Epetra_Vector> RestrictVector(Epetra_Vector const &vector,
                                               std::vector<int> const &indices);


    //! a general function for gathering sparse matrices. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    Teuchos::RCP<Epetra_CrsMatrix> Gather(const Epetra_CrsMatrix& mat, int root);

    //! a general function for gathering vectors. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    Teuchos::RCP<Epetra_MultiVector> Gather(const Epetra_MultiVector& vec, int root);

    //! transform a "solve" or "standard" into a replicated "gather" map
    //! The new map will have its indices sorted in ascending order.
    Teuchos::RCP<Epetra_BlockMap> Gather(const Epetra_BlockMap& map, int root);

    //! transform a "solve" or "standard" into a replicated "gather" map
    //! The new map will have its indices sorted in ascending order.
    Teuchos::RCP<Epetra_Map> Gather(const Epetra_Map& map, int root);

    //! convert a 'gathered' vector into a distributed vector
    Teuchos::RCP<Epetra_MultiVector> Scatter(const Epetra_MultiVector& vec,
                                             const Epetra_BlockMap& distmap);

    //! transform a "solve" into a replicated "col" map, i.e.
    //! the resulting map will have all gid's of the original one
    //! on every process. The new map will have its indices
    //! sorted in ascending order if you choose reorder=true,
    //! otherwise the ordering is retained as it is.
    Teuchos::RCP<Epetra_BlockMap> AllGather(const Epetra_BlockMap& map,
                                            bool reorder=true);

    //! transform a "solve" into a replicated "col" map, i.e.
    //! the resulting map will have all gid's of the original one
    //! on every process. The new map will have its indices
    //! sorted in ascending order if you choose reorder=true,
    //! otherwise the ordering is retained as it is.
    Teuchos::RCP<Epetra_Map> AllGather(const Epetra_Map& map, bool reorder=true);

    //! a general function for gathering vectors. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    Teuchos::RCP<Epetra_MultiVector> AllGather(const Epetra_MultiVector& vec);

    //! a general function for gathering vectors. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    Teuchos::RCP<Epetra_IntVector> AllGather(const Epetra_IntVector& vec);

    //! compute matrix-matrix product C=A*B (implemented using EpetraExt)
    Teuchos::RCP<Epetra_CrsMatrix> MatrixProduct(bool transA, const Epetra_CrsMatrix& A,
                                                 bool transB, const Epetra_CrsMatrix& B,
                                                 bool useColMap = false);

    //! compute triple matrix product A*B*C or the like
    Teuchos::RCP<Epetra_CrsMatrix> TripleProduct(bool transA, const Epetra_CrsMatrix& A,
                                                 bool transB, const Epetra_CrsMatrix& B,
                                                 bool transC, const Epetra_CrsMatrix& C);

    //! extracts those nodes from a column map that are actually present in a
    //! Matrix' Graph
    Teuchos::RCP<Epetra_Map> CompressColMap(const Epetra_CrsMatrix& A);

    //! generalized ReplaceColMap that only requires the new col map to be a
    //! subset of the old one (instead of being 'PointSameAs' as in the Trilinos
    //! function). The map should contain all column indices actually appearing,
    //! such a map can be created by 'CompressColMap'.
    Teuchos::RCP<Epetra_CrsMatrix> ReplaceColMap(Teuchos::RCP<Epetra_CrsMatrix> A,
                                                 const Epetra_Map& newmap);
    //! Workaround for the buggy 'ReplaceRowMap' routine in Trilinos 7.0.
    /*! It is much slower as it copies the whole matrix and should only
        be used until they fix the builtin one
    */
    Teuchos::RCP<Epetra_CrsMatrix> ReplaceRowMap(Teuchos::RCP<Epetra_CrsMatrix> A,
                                                 const Epetra_Map& newmap);

    //! create an exact copy of a matrix removing the column map.
    //! This means that row- and column map have to be 'compatible'
    //! in the sense that the ColMap is a subset of the RowMap.
    //! It seems to be required in order to use Ifpack in some cases.
    Teuchos::RCP<Epetra_CrsMatrix> RemoveColMap(Teuchos::RCP<Epetra_CrsMatrix> A);

    Teuchos::RCP<Epetra_CrsMatrix> RebuildMatrix(Teuchos::RCP<Epetra_CrsMatrix> A);

    //! simultaneously replace row and column map (see comment for previous function)
    //! This is a special purpose function. The newcolmap must be a subset of the current
    //! colmap, i.e. you cannot really change the indexing scheme for the columns.
    Teuchos::RCP<Epetra_CrsMatrix> ReplaceBothMaps(Teuchos::RCP<Epetra_CrsMatrix> A,
                                                   const Epetra_Map& newmap,
                                                   const Epetra_Map& newcolmap);

}
#endif
