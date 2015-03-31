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

// A namespace with the necessary utils from hymls and trios
// It should at least contain:
//
//   HYMLS::Tools::SplitBox
//   HYMLS::Tools::ind2sub
//   
//   HYMLS::MatrixUtils::CreateMap
//   HYMLS::MatrixUtils::Gather
//   HYMLS::MatrixUtils::Scatter

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
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
}
#endif
