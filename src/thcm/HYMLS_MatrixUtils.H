/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include "HYMLS_Macros.H"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "AnasaziTypes.hpp"

// this is useful for comparing with MATLAB:
//#define REINDEX_BY_DEFAULT true

// for comparing internally (e.g. with different number of
// processors) it is better not to re-index when printing 
// vectors or matrices:
#ifdef MATLAB_COMPATIBILITY_MODE
#define REINDEX_BY_DEFAULT true
#else
#define REINDEX_BY_DEFAULT false
#endif
class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_Operator;
class Epetra_Map;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_Import;
class Ifpack_Container;

namespace HYMLS {

//! utility class for the block-ILU preconditioner

/*! a purely static class that offers auxiliary routines
     for splitting and reordering matrices/vectors and
     related operations.
     */
class MatrixUtils
  {
  
  public:
  
  typedef enum {
  Absolute,   // drop aij if abs(aij) <= tol
  AbsZeroDiag,// like Absolute, but set diagonal entries with abs(aii)<tol to 0.0 rather 
              // than discarding them
  Relative, // drop aij if abs(aij) <= tol*max(|aii|,|ajj|). 
            // This behavior should prevent asymmetric dropping in F-matrices,
            // diagonal entries are never dropped (unless you set tol>=1).
  RelDropDiag, // like Relative, but use absolute dropping criterion on diagonal,
               // e.g. delete aii if |aii|<=tol
  RelZeroDiag // like RelDropDiag, but put physical 0.0 in aii instead of deleting it.
  } DropType;

  typedef enum {
  MATRIXMARKET, // use EpetraExt ToMatrixMarketFile functions (fails for overlapping objects)
  GATHER,       // collect to root and dump in ascii file (may require lots of memory on one proc)
  HDF5          // not implemented, probably the best thing to do as it is portable
  } PrintMethod;
  
  typedef Anasazi::Eigensolution<double,Epetra_MultiVector> Eigensolution;

    //! create a map for a 3D problem with multiple degrees of freedom per node
    
    /*! creates a distributed non-overlapping map to be used for constructing/solving
        linear problems with one or more dof/node. This is similar to the Galeri map
        Cartesian3D, except that we allow several dof/node.
        For 2D problems simply set nz=1.
        
        This function can also be used if the number of grid cells nx*ny*nz
        is less than the number of processes, in which case any proc with PID
        >= nx*ny*nz gets 0 elements in the map.
    
    It is possible to artifically restrict the number of procs that get any nodes to
    maxNumActiveProcs by setting the flag to something >0. The default -1 means all in comm 
    get nodes (if there are sufficient nodes).
    */    
    static Teuchos::RCP<Epetra_Map> CreateMap(int nx, int ny, int nz, 
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
    static Teuchos::RCP<Epetra_Map> CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
                             int I0, int I1, int J0, int J1, int K0, int K1,
                             int dof,
                             const Epetra_Comm& comm);
  
    //! extract a map with 1 variable from a map with dof variables per node. 'var'
    //! is the variable to be extracted, i.e. 0, 1,...dof-1 
    static Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, int dof, int var);

    //! extract a map with all indices of the original one between min_gid and max_gid 
    //! (including the bounds)
    static Teuchos::RCP<Epetra_Map> ExtractRange(const Epetra_Map& map, int min_gid, int max_gid);

    //! extract a map for two variables from a map with dof variables. 'var'
    //! are the variables to be extracted, i.e. {0,1}, {4,5} etc.
    static Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, int dof, const int var[2]);

    //! extract a map for nvars variables from a map with dof variables. 'var'
    //! are the variables to be extracted.
    static Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, int dof, const int *var, int nvars);
    
    //! given a map and an array indicating wether each node of the map is to be
    //! discarded (true) or not (false), this function creates a new map with the
    //! discarded entries removed.
    static Teuchos::RCP<Epetra_Map> CreateSubMap(const Epetra_Map& map, const bool* discard);    
                        

    //! extracts those nodes from a column map that are actually present in a
    //! Matrix' Graph
    static Teuchos::RCP<Epetra_Map> CompressColMap(const Epetra_CrsMatrix& A);
    
    //! create an optimal column map for extracting A(rowMap,colMap), given a distributed
    //! column map which has entries owned by other procs that we need for the column map.
    static Teuchos::RCP<Epetra_Map> CreateColMap(const Epetra_CrsMatrix& A, 
        const Epetra_Map& newRows, const Epetra_Map& newColumns);

    //! a general function for gathering sparse matrices. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    static Teuchos::RCP<Epetra_CrsMatrix> Gather(const Epetra_CrsMatrix& mat, int root);

    //! a general function for gathering vectors. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    static Teuchos::RCP<Epetra_MultiVector> Gather(const Epetra_MultiVector& vec, int root);

    //! a general function for gathering vectors. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    static Teuchos::RCP<Epetra_IntVector> Gather(const Epetra_IntVector& vec, int root);

    //! transform a "solve" or "standard" into a replicated "gather" map
    //! The new map will have its indices sorted in ascending order.    
    static Teuchos::RCP<Epetra_BlockMap> Gather(const Epetra_BlockMap& map, int root);

    //! transform a "solve" or "standard" into a replicated "gather" map
    //! The new map will have its indices sorted in ascending order.    
    static Teuchos::RCP<Epetra_Map> Gather(const Epetra_Map& map, int root);

    //! transform a "solve" into a replicated "col" map, i.e.     
    //! the resulting map will have all gid's of the original one 
    //! on every process. The new map will have its indices       
    //! sorted in ascending order if you choose reorder=true,     
    //! otherwise the ordering is retained as it is.              
    static Teuchos::RCP<Epetra_BlockMap> AllGather(const Epetra_BlockMap& map, bool reorder=true);

    //! transform a "solve" into a replicated "col" map, i.e.     
    //! the resulting map will have all gid's of the original one 
    //! on every process. The new map will have its indices       
    //! sorted in ascending order if you choose reorder=true,     
    //! otherwise the ordering is retained as it is.              
    static Teuchos::RCP<Epetra_Map> AllGather(const Epetra_Map& map, bool reorder=true);

    //! a general function for gathering vectors. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    static Teuchos::RCP<Epetra_MultiVector> AllGather(const Epetra_MultiVector& vec);

    //! a general function for gathering vectors. It is fairly slow
    //! as it rebuilds the required "GatherMap" every time.
    static Teuchos::RCP<Epetra_IntVector> AllGather(const Epetra_IntVector& vec);
    
    //! convert a 'gathered' vector into a distributed vector
    static Teuchos::RCP<Epetra_MultiVector> Scatter(const Epetra_MultiVector& vec, 
                                               const Epetra_BlockMap& distmap);

    //! Workaround for the buggy 'ReplaceRowMap' routine in Trilinos 7.0.
    
    /*! It is much slower as it copies the whole matrix and should only 
        be used until they fix the builtin one
    */
    static Teuchos::RCP<Epetra_CrsMatrix> ReplaceRowMap(Teuchos::RCP<Epetra_CrsMatrix> A, const Epetra_Map& newmap);

    //! generalized ReplaceColMap that only requires the new col map to be a 
    //! subset of the old one (instead of being 'PointSameAs' as in the Trilinos
    //! function). The map should contain all column indices actually appearing,
    //! such a map can be created by 'CompressColMap'.
    static Teuchos::RCP<Epetra_CrsMatrix> ReplaceColMap(Teuchos::RCP<Epetra_CrsMatrix> A, const Epetra_Map& newmap);

  //! create an exact copy of a matrix removing the column map.
  //! This means that row- and column map have to be 'compatible'
  //! in the sense that the ColMap is a subset of the RowMap.
  //! It seems to be required in order to use Ifpack in some cases.
  static Teuchos::RCP<Epetra_CrsMatrix> RemoveColMap(Teuchos::RCP<Epetra_CrsMatrix> A);

    //! simultaneously replace row and column map (see comment for previous function)
    //! This is a special purpose function. The newcolmap must be a subset of the current
    //! colmap, i.e. you cannot really change the indexing scheme for the columns.
    static Teuchos::RCP<Epetra_CrsMatrix> ReplaceBothMaps(Teuchos::RCP<Epetra_CrsMatrix> A, const Epetra_Map& newmap, 
                                            const Epetra_Map& newcolmap);
    
    //! local upper or lower tri-solve. This can be used to gain insight into the process,
    //! but in principle the Trilinos version should work fine.
    static void TriSolve(const Epetra_CrsMatrix& A,const Epetra_Vector& b, Epetra_Vector& x);
        
    //! replace matrix by identity matrix
    
    /*! All matrix entries are set to 0 and the diagonal is set to 1.
    this is used for testing the block ILU precondioner step-wise and
    is not at all foolproof (it may fail when there is no diagonal   
    in the original matrix, for instance)
    */
    static void Identity(Teuchos::RCP<Epetra_CrsMatrix> A);
    

    //! dump CrsMatrix to file (the matrix is gathered so it can be easily
    //! read into MATLAB etc., but this is only meant for debugging etc.)
    static void Dump(const Epetra_CrsMatrix& A, const std::string& filename,
           bool reindex=REINDEX_BY_DEFAULT,
           PrintMethod how=MATRIXMARKET);
    
    //! dump CrsMatrix in binary file (HDF5 format), which is parallel,   
    //! space efficient and preserves the data distribution. If the code  
    //! is compiled without -DHAVE_XDMF, this will cause an error.        
    //! The groupname is used to locate the stored matrix within the HDF- 
    //! file.                                                             
    //! if new_file==true, a new HDF file is created, deleting any old    
    //! file with the same name. Otherwise, an old file is opened, which  
    //! has to exist or an exception will be thrown.                      
    static void DumpHDF(const Epetra_CrsMatrix& A, const std::string& filename, 
        const std::string& groupname, bool new_file=true);

    //! dump Vector to file (the vector is gathered so it can be easily
    //! read into MATLAB etc., but this is only meant for debugging etc.)
    static void Dump(const Epetra_MultiVector& x, const std::string& filename, bool 
        reindex=REINDEX_BY_DEFAULT,
        PrintMethod how=GATHER);

    //! dump Vector to file (the vector is gathered so it can be easily
    //! read into MATLAB etc., but this is only meant for debugging etc.)
    static void Dump(const Epetra_IntVector& x, const std::string& filename);

    //! dump Map to file (the map is gathered so it can be easily
    //! read into MATLAB etc., but this is only meant for debugging etc.)
    static void Dump(const Epetra_Map& M, const std::string& filename,
        PrintMethod how=MATRIXMARKET);
        
    //! dump vector in binary (HDF5) format (see comment on DumpMatrixHDF)
    static void DumpHDF(const Epetra_MultiVector& x, const std::string& filename, 
        const std::string& groupname, bool new_file=true);
    
    
    //! Print an Epetra_RowMatrix to a stream
    static void PrintRowMatrix(const Epetra_RowMatrix&,std::ostream&);
    
    //! MatrixMarket output of MultiVector. The vectors are written in order
    //! of ascending GID, but without storing the GIDs. This is different from
    //! just calling the EpetraExt function as that one stores the entries    
    //! processor-by-processor.
    static int mmwrite(std::string filename, const Epetra_MultiVector& vec);

    //! MatrixMarket input of MultiVector.
    static int mmread(std::string filename, Epetra_MultiVector& vec);
    
    //! compute triple matrix product A*B*C or the like
    static Teuchos::RCP<Epetra_CrsMatrix> TripleProduct(bool transA, const Epetra_CrsMatrix& A,
                                 bool transB, const Epetra_CrsMatrix& B,
                                 bool transC, const Epetra_CrsMatrix& C);

    //! compute triple matrix product M=A*B*C or the like if the pattern of M is known
    static void TripleProduct(Teuchos::RCP<Epetra_CrsMatrix> M, bool transA, const Epetra_CrsMatrix& A,
                                 bool transB, const Epetra_CrsMatrix& B,
                                 bool transC, const Epetra_CrsMatrix& C);



    //! compute matrix-matrix product C=A*B (implemented using EpetraExt)
    static Teuchos::RCP<Epetra_CrsMatrix> MatrixProduct(bool transA, const Epetra_CrsMatrix& A,
                                 bool transB, const Epetra_CrsMatrix& B);

    //! compute matrix-matrix product C=A*B (implemented using EpetraExt) if the patern of C 
    //! is already known (C is Filled)
    static void MatrixProduct(Teuchos::RCP<Epetra_CrsMatrix> M, bool transA, const Epetra_CrsMatrix& A,
                                 bool transB, const Epetra_CrsMatrix& B);

    //! compute a number of eigenvalues and eigenvectors of the sparse
    //! matrix pencil [A,B]. The most dominant eigs are found, so if  
    //! you want e.g. those of smallest magnitude, you have to supply 
    //! thie inverse operator yourself. The solution is returned as an
    //! Anasazi::Eigensolution, which contains eigenvalues, -vectors  
    //! and an orthogonal basis for the found eigenspace.             
    //!                                                               
    //! The method used is the Anasazi Block-Krylov-Schur algorithm   
    //! (Arnoldi, no symmetry is taken into account).                 
    //!                                                               
    //! B may be Teuchos::null, in which case the standard eigenvalue 
    //! problem Ax=lambda x is solved.                                
    //!                                                               
    //! The eigenvectors returned are the right eigenvectors of [A,B].
    //!                                                               
    static Teuchos::RCP<Eigensolution> Eigs(
                Teuchos::RCP<const Epetra_Operator> A, 
                Teuchos::RCP<const Epetra_Operator> B, 
                int howMany=6,
                double tol=1.0e-6);
                                 
    //! reads sparse matrix from THCM's beg/co/jco/info format
    //! if the prefix is, for instance, "A", then the matrix  
    //! is read from A.info, A.beg etc.
    static Teuchos::RCP<Epetra_CrsMatrix> ReadThcmMatrix(std::string prefix, const Epetra_Comm& comm, 
           const Epetra_Map& rowmap, const Epetra_Map& colmap,
           const Epetra_Map* rangemap=NULL,const Epetra_Map* domainmap=NULL);
                       
    //! this is for creating consistent random vectors. Normally, if you change the number 
    //! of procs and use the same seed in a v.Random() call, the vector is different. 
    //! Actually I think the random vectors generated that way are highly correlated or
    //! even identical on every processor, so we just collect the vector, randomize it and
    //! distribute it again (note that this function is really inefficient, it is just
    //! intended for creating initial vectors for reproducible results)
    static int Random(Epetra_MultiVector& v, int seed=-1);
    
    //! drop elements <= a threshold from a sparse matrix. By default the 
    //! threshold is the machine epsilon so that numerical zeros are dropped.
    //! The threshold is scaled in a certain way, specified by the t parameter.
    //! For details, see the DropType enum.
    static Teuchos::RCP<Epetra_CrsMatrix> DropByValue(Teuchos::RCP<const Epetra_CrsMatrix> A, 
        double threshold=HYMLS_SMALL_ENTRY, DropType t=RelZeroDiag);
    
    //! replace one row and column by a Dirichlet condition (0 everywhere, 1 on diagonal).
    //! This function assumes that the pattern of the matrix is symmetric!
    static int PutDirichlet(Epetra_CrsMatrix& A, int gid);
    
    
    //!                                                                         
    //! computes a fill-reducing ordering for a sparse matrix K.                
    //!                                                                         
    //! - For general sparse matrices with non-zero diagonal, metis is used to  
    //!   create the permutation (via Isorropia and Zoltan), and                
    //!   row_perm=col_perm.                                                    
    //! - for serial F-matrices, an orderinbased on the pattern of A+BB' is     
    //!   computed and the P-nodes are inserted according to Fred and Arie's    
    //!   algorithm. The matrix K(row_perm,col_perm) will have no zero diagonal 
    //!   entries.                                                              
    //! - otherwise currently a warning is issued and an identity is returned   
    //!   as ordering (TODO: we need the F-mat ordering also for the reduced SC)
    //!                                                                         
    //! if dummy is true, no actual fill-reducing ordering is computed. In case 
    //! of an F-matrix, the p-nodes are inserted correctly into the natural or- 
    //! dering. dummy defaults to false.                                        
    //!                                                                         
    static int FillReducingOrdering(const Epetra_CrsMatrix& K,
                                Teuchos::Array<int>& row_perm,
                                Teuchos::Array<int>& col_perm,
                                bool dummy=false);

    //! computes the AMD ordering of a serial input matrix using AMD from SuiteSparse.
    static int AMD(const Epetra_CrsGraph& A,
                        Teuchos::Array<int>& perm);

    //! creates a map of the form (i0:i1,j0:j1,k0:k1),
    //! column-major access (i is fastest).           
    //! Small letters indicate the bounds of the      
    //! subdomain, capital letters the bounds of the  
    //! global domain. The indices of subdomains are  
    //! allowed to overlap                          */
    static Teuchos::RCP<Epetra_Map> CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
                             int I0, int I1, int J0, int J1, int K0, int K1,
                             const Epetra_Comm& comm);
    

    //! sort an integer and a double array consistently so the integer array
    //! is in ascending order
    static int SortMatrixRow(int* inds, double* vals, int len);

    //! extract a local part of a matrix. This can not easily be done by Import
    //! objects because they tend to do collective communication. A_loc should 
    //! be !DistributedGlobal(), and A should be Filled (if it is Crs).        
    //! We do not call FillComplete in this function.
    static int ExtractLocalBlock(const Epetra_RowMatrix& A, Epetra_CrsMatrix& A_loc);

  private:
  
    //! returns a string describing the class
    static std::string Label() {return "MatrixUtils";}
  
    //! read array from file
    static void read_fortran_array(int n, int* array, std::string filename);
    
    //! read array from file
    static void read_fortran_array(int n, double* array, std::string filename);

  };

}

#endif // MATRIX_UTILS_h

