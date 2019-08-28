/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef OCEANGRID_H
#define OCEANGRID_H

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

#include "THCMdefs.H"

#include "TRIOS_Domain.H"

#include "Teuchos_Array.hpp"

typedef int MaskType;

//! Grid representation of model vectors

/*! representation of a model vector containing u,v,w,p,T,S
  as a set of 3D fields. This class assumes the grid looks
  like the one used in THCM, i.e.

  ^           +---w---+        +______uv  ^
  |           |       |        |       |  |
  z           |  pTS  uv       | wpTS  | phi
  |           |       |        |       |  |
  theta->     +-------+        +-------+ -theta->

  (front view)      (top view)

  i: theta
  j: phi
  k: z

  the index ranges for the fields are:

  u,v:   i=0:n,   j=0:m,   k=0:l+1
  w:     i=0:n+1, j=0:m+1, k=0:l
  p,T,S: i=0:n+1, j=0:m+1, k=0:l+1
*/
class OceanGrid
{
public:
    //!
    OceanGrid(Teuchos::RCP<TRIOS::Domain> dom);

    //!
    virtual ~OceanGrid();

    //! put data from vector in grid format.
    //! input should be based on the 'Assembly' map
    void ImportData(const Epetra_Vector& input);

    //! import a land mask. Should be an overlapping object
    //! already, i.e. (0:n+1,0:m+1,0:l+1)
    void ImportLandMask(const Epetra_IntVector& mask);

    //! put data from grid in a vector.
    //! input should be based on the 'Assembly' map
    void ExportData(Epetra_Vector& output);

    // output to a stream (you can use operator << instead)
    std::ostream& print(std::ostream& os) const;


    //! access U-data
    inline double& u(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_UV(i,j,k);
#endif
            return U_[index_UV(i,j,k)];
        }

    //! access U-data (const)
    inline const double& u(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_UV(i,j,k);
#endif
            return U_[index_UV(i,j,k)];
        }

    //! access V-data
    inline double& v(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_UV(i,j,k);
#endif
            return V_[index_UV(i,j,k)];
        }

    //! access V-data (const)
    inline const double& v(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_UV(i,j,k);
#endif
            return V_[index_UV(i,j,k)];
        }

    //! access U-data
    inline double& w(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_W(i,j,k);
#endif
            return W_[index_WPTS(i,j,k)];
        }

    //! access W-data (const)
    inline const double& w(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_W(i,j,k);
#endif
            return W_[index_WPTS(i,j,k)];
        }

    //! access U-data
    inline double& p(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return P_[index_WPTS(i,j,k)];
        }

    //! access P-data (const)
    inline const double& p(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return P_[index_WPTS(i,j,k)];
        }

    //! access T-data
    inline double& T(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return T_[index_WPTS(i,j,k)];
        }

    //! access T-data (const)
    inline const double& T(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return T_[index_WPTS(i,j,k)];
        }

    //! access S-data
    inline double& S(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return S_[index_WPTS(i,j,k)];
        }

    //! access S-data (const)
    inline const double& S(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return S_[index_WPTS(i,j,k)];
        }

    //! access density-data
    inline double& rho(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return S_[index_WPTS(i,j,k)];
        }

    //! access density-data (const)
    inline const double& rho(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return Rho_[index_WPTS(i,j,k)];
        }

    //! access land mask (range is 0:m+1,0:n+1,0:l+1)
    MaskType landm(int i,int j,int k)
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return LandMask_[index_WPTS(i,j,k)];
        }

    //! access land mask (range is 0:m+1,0:n+1,0:l+1)
    MaskType landm(int i,int j,int k) const
        {
#ifdef TESTING
            bounds_check_PTS(i,j,k);
#endif
            return LandMask_[index_WPTS(i,j,k)];
        }

    //! get domain-decomposition object
    Teuchos::RCP<TRIOS::Domain> getDomain() {return domain;}

    //! returns the maximum of the meridional streamfunction Psi_M.
    //! If necessary, it is recomputed.
    inline double psimMax(void)
        {
            if (recompute_PsiM_) recomputePsiM();
            return PsimMax_;
        }

    //! returns the minimum of the meridional streamfunction Psi_M.
    //! If necessary, it is recomputed.
    inline double psimMin(void)
        {
            if (recompute_PsiM_) recomputePsiM();
            return PsimMin_;
        }

    //! returns the maximum of the barotropic streamfunction Psi_B.
    //! If necessary, it is recomputed.
    inline double psibMax(void)
        {
            if (recompute_PsiB_) recomputePsiB();
            return PsibMax_;
        }

    //! returns the minimum of the barotropic streamfunction Psi_M.
    //! If necessary, it is recomputed.
    inline double psibMin(void)
        {
            if (recompute_PsiB_) recomputePsiB();
            return PsibMin_;
        }

    //! returns the maximum of the horizontal velocity u
    //! (absolute value) evaluated at the cell center.
    //! If necessary, it is recomputed.
    inline double uMax(void)
        {
            if (recompute_MaxVel_) recomputeMaxVel();
            return MaxU;
        }

    //! returns the maximum of the horizontal velocity v
    //! (absolute value) evaluated at the cell center.
    //! If necessary, it is recomputed.
    inline double vMax(void)
        {
            if (recompute_MaxVel_) recomputeMaxVel();
            return MaxV;
        }

    //! returns the maximum of the vertical velocity w
    //! (absolute value) evaluated at the cell center.
    //! If necessary, it is recomputed.
    inline double wMax(void)
        {
            if (recompute_MaxVel_) recomputeMaxVel();
            return MaxW;
        }

    double cflCond(void)
        {
            return std::min(std::min(dx/uMax(),dy/vMax()),dz/wMax());
        }

    void GetPointers(double **u, double **v, double **w,
                     double **p, double **t, double **s)
        {
            *u = U_; *v = V_; *w = W_;
            *p = P_; *t = T_; *s = S_;
        }
    //! \name get coordinates of cell centers
    //@{
    //! get x coordinates
    inline Teuchos::RCP<Epetra_Vector> GetXc(){return xc_;}
    //! get y coordinates
    inline Teuchos::RCP<Epetra_Vector> GetYc(){return yc_;}
    //! get z coordinates
    inline Teuchos::RCP<Epetra_Vector> GetZc(){return zc_;}
    //@}

    //! \name get coordinates of grid nodes
    //@{
    //! get x coordinates
    inline Teuchos::RCP<Epetra_Vector> GetXu(){return xu_;}
    //! get y coordinates
    inline Teuchos::RCP<Epetra_Vector> GetYv(){return yv_;}
    //! get z coordinates
    inline Teuchos::RCP<Epetra_Vector> GetZw(){return zw_;}
    //@}

protected:

    //! object describing the domain-decomposition
    Teuchos::RCP<TRIOS::Domain> domain;

    //! vector for importing ghost-values into the subdomain
    Teuchos::RCP<Epetra_Vector> importVector;

    //! these maps are for one-dimensional quantities, i.e. the
    //! arrays containing x/y/z coordinates
    Teuchos::RCP<Epetra_Map> xMap,yMap,zMap;

    //! these maps are 1D 0-based maps (they have one additional 'cell' i/j/k=0)
    Teuchos::RCP<Epetra_Map> xMap0,yMap0,zMap0;

    //! vectors containing the coordinates of the cell-centers
    Teuchos::RCP<Epetra_Vector> xc_,yc_,zc_;

    //! vectors containing the coordinates of the nodes (cell-corners)
    Teuchos::RCP<Epetra_Vector> xu_,yv_,zw_;

    //! communicators only in x- and y-direction
    Teuchos::RCP<Epetra_Comm> xComm, yComm;

    //! data arrays for the six unknowns
    double *U_, *V_, *W_, *P_, *T_, *S_;

    //! grid spacing
    double dx,dy,dz;

    //! array containing the meridional streamfunction
    //! (this array is not yet accessible, only PsiMax
    //! can be requested whereupon PsiM is computed)
    double *PsiM_, *PsiB_, *Rho_;

    //! landm array (copy of the one in Fortran)
    MaskType *LandMask_;


    //! global min/maximum of the meridional streamfunction
    double PsimMin_,PsimMax_;

    //! global min/maximum of the barotropic streamfunction
    double PsibMin_,PsibMax_;

    //! maximum (absolute) of velocity components at cell-centers
    double MaxU,MaxV,MaxW;

    //! local array sizes, note that this includes ghost cells
    int n, m, l;

    // after a new import: which quatities are there, which have to be recomputed?
    bool recompute_PsiM_, recompute_PsiB_, recompute_MaxVel_;


private:

    inline void bounds_check_UV(int i,int j,int k) const
        {
            if ((i<0)||(j<0)||(k<0)||(i>n)||(j>m)||(k>l+1))
            {
                INFO("index out of bounds, u/v array, index: ("<<i<<", "<<j<<", "<<k<<")");
                INFO("                         array bounds: (0:"<<n<<", 0:"<<m<<", 0:"<<l+1<<")");
                ERROR("Index out of bounds!",__FILE__,__LINE__);
            }
        }

    inline void bounds_check_W(int i,int j,int k) const
        {
            if ((i<0)||(j<0)||(k<0)||(i>n+1)||(j>m+1)||(k>l))
            {
                INFO("index out of bounds, w array, index: ("<<i<<", "<<j<<", "<<k<<")");
                INFO("                       array bounds: (0:"<<n+1<<", 0:"<<m+1<<", 0:"<<l<<")");
                ERROR("Index out of bounds!",__FILE__,__LINE__);
            }
        }

    inline void bounds_check_PTS(int i,int j,int k) const
        {
            if ((i<0)||(j<0)||(k<0)||(i>n+1)||(j>m+1)||(k>l+1))
            {
                INFO("index out of bounds, p/T/S array, index: ("<<i<<", "<<j<<", "<<k<<")");
                INFO("                       array bounds: (0:"<<n+1<<", 0:"<<m+1<<", 0:"<<l+1<<")");
                ERROR("Index out of bounds!",__FILE__,__LINE__);
            }
        }

    inline void bounds_check_PsiM(int j,int k) const
        {
            if ((j<0)||(k<0)||(j>m)||(k>l))
            {
                INFO("index out of bounds, psiM array, index: ("<<j<<", "<<k<<")");
                INFO("                       array bounds: (0:"<<m<<", 0:"<<l<<")");
                ERROR("Index out of bounds!",__FILE__,__LINE__);
            }
        }

    inline void bounds_check_PsiB(int i,int j) const
        {
            if ((i<0)||(j<0)||(i>n)||(j>m))
            {
                INFO("index out of bounds, psiB array, index: ("<<i<<", "<<j<<")");
                INFO("                       array bounds: (0:"<<n<<", 0:"<<m<<")");
                ERROR("Index out of bounds!",__FILE__,__LINE__);
            }
        }

    inline int index_UV(int i, int j, int k) const
        {
            return ((m+1)*k + j)*(n+1) + i;
        }

    inline int index_WPTS(int i, int j, int k) const
        {
            return ((m+2)*k + j)*(n+2) + i;
        }


    inline int index_PsiM(int j, int k) const
        {
            return (m+1)*k + j;
        }

    inline int index_PsiB(int i, int j) const
        {
            return (n+1)*j + i;
        }


    //! access meridional streamfunction (may become public at some point)
    inline const double& psiM(int j,int k) const
        {
#ifdef TESTING
            bounds_check_PsiM(j,k);
#endif
            return PsiM_[index_PsiM(j,k)];
        }

    //! access meridional streamfunction (may become public at some point)
    inline double& psiM(int j,int k)
        {
#ifdef TESTING
            bounds_check_PsiM(j,k);
#endif
            return PsiM_[index_PsiM(j,k)];
        }



    //! access barotropic streamfunction (may become public at some point)
    inline const double& psiB(int i,int j) const
        {
#ifdef TESTING
            bounds_check_PsiB(i,j);
#endif
            return PsiB_[index_PsiB(i,j)];
        }

    //! access meridional streamfunction (may become public at some point)
    inline double& psiB(int i,int j)
        {
#ifdef TESTING
            bounds_check_PsiB(i,j);
#endif
            return PsiB_[index_PsiB(i,j)];
        }




    //! memory allocation helper function

    /*! all arrays are handled as 1D arrays internally,
      giving 3D array bounds makes the function independent
      of actual indexing schemes etc.
    */
    double* allocateGrid(int imin, int imax,
                         int jmin, int jmax,
                         int kmin, int kmax);

    //! recompute the meridional streamfunction PsiM
    void recomputePsiM(void);

    //! recompute the meridional streamfunction PsiM
    void recomputePsiB(void);

    //! recompute the maximum velocities
    void recomputeMaxVel(void);

};

//! output operator
std::ostream& operator<<(std::ostream& os, const OceanGrid& G);


#endif
