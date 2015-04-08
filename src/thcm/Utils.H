/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
/**********************************************************************
 * Put together by T.E. Mulder, Utrecht University 2014/15            *
 * contact: t.e.mulder@uu.nl                                          *
 **********************************************************************/
#ifndef UTILS_H
#define UTILS_H

// A namespace with utils from hymls and trios that I need

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

// for sorting indices
#include <algorithm>

#include "GlobalDefinitions.H"

namespace Utils
{
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
												 bool transB, const Epetra_CrsMatrix& B);

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

	//! simultaneously replace row and column map (see comment for previous function)
    //! This is a special purpose function. The newcolmap must be a subset of the current
    //! colmap, i.e. you cannot really change the indexing scheme for the columns.
    Teuchos::RCP<Epetra_CrsMatrix> ReplaceBothMaps(Teuchos::RCP<Epetra_CrsMatrix> A,
												   const Epetra_Map& newmap,
												   const Epetra_Map& newcolmap);

}
#endif
