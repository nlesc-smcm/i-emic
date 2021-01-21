/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#ifndef TRIOS_DOMAIN_H
#define TRIOS_DOMAIN_H

#include <Teuchos_ArrayRCPDecl.hpp>
#include <Teuchos_RCP.hpp>


#include <functional>
#include <vector>

class Epetra_Comm;
class Epetra_Map;
class Epetra_Vector;
class Epetra_CrsMatrix;
class Epetra_Import;

namespace TRIOS {
    //! (sub)domain geometry
    struct Grid
    {
        int n_, m_, l_;
        double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_, hdim_;

        Teuchos::ArrayRCP<double> x_, y_, z_, xu_, yv_, zw_;

        Grid();
        Grid(double qz, int N, int M, int L, double xmin, double xmax,
             double ymin, double ymax, double hdim);

        Grid SubGrid(int Noff, int Moff, int nloc, int mloc, int lloc) const;

        double getXposEdge(int x) const;
        double getXposCenter(int x) const;
        double getYposEdge(int y) const;
        double getYposCenter(int y) const;
        double getZposEdge(int z) const;
        double getZposCenter(int z) const;

      private:
        // vertical grid stretching function
        double fz(double);
        double dfdz(double);

        //! vertical grid stretching
        double qz_;

        Teuchos::ArrayRCP<double> depth_, cumulativeDepth_;
    };

    //! domain decomposition class

    /*! This is a class to define the geometry of a rectangular
      (ocean) domain.
      The main purpose of this class is to decompose the domain
      and create Epetra_Map objects (function DecompXD).
      There are two maps, one including ghost-nodes (Assembly Map),
      whereas the other (Solve Map) is a 1-to-1 mapping.
      Corresponding Export/Import objects are constructed as well.
    */
    class Domain
    {

    public:
        //! Constructor
        /*!
          input: dimensions of the global box:
          \verbatim
          - N: east-west
          - M: north-south
          - L: z-direction
          \endverbatim
        */
        Domain(int N, int M, int L, int dof,
               double xmin, double xmax, double ymin, double ymax,
               bool periodic, double hdim, double qz,
               Teuchos::RCP<Epetra_Comm> comm, int aux = 0);

        //! Destructor
        virtual ~Domain();

        //! decompose the domain for a 2D processor array.
        /*! The xy-directions are split up.

          The domain is split up in a 2D column-major fashion, i.e.
          (we call it column-major because it comes from the Fortran THCM impl.)

          \verbatim

          _______ ^
          |P1|P2| |
          +--+--+ m
          |P3|P4| |
          +--+--+ v
          <--n-- >

          \endverbatim

          Ghost-nodes are added between subdomains, i.e.
          grid point (i,mloc,k) on P2 ^= (i,1,k) on P1 etc.
          Two maps are created, one including ghost-nodes (the assembly map),
          the other not including ghost-nodes (the solve map).
        */
        void Decomp2D();

        //! return number of grid cells in the subdomain
        //! in east-west direction (including ghost-nodes)
        inline int LocalN() const {return nloc;}

        //! return number of grid points in the subdomain
        //! in north-south direction (including ghost-nodes)
        inline int LocalM() const {return mloc;}

        //! return number of grid points in the subdomain
        //! in z direction (including ghost-nodes)
        inline int LocalL() const {return lloc;}

        //! return number of grid points in the entire domain
        //! in east-west direction
        inline int GlobalN() const {return gridGlb_.n_;}

        //! return number of grid points in the entire domain
        //! in north-south direction
        inline int GlobalM() const {return gridGlb_.m_;}

        //! return number of grid points in the entire domain
        //! in z direction
        inline int GlobalL() const {return gridGlb_.l_;}

        //! return dof
        inline int Dof() const {return dof_;}

        //! return aux
        inline int Aux() const {return aux_;}

        //! get first global k index (z-direction) on the subdomain
        //! (including ghost-nodes)
        inline int FirstK(){return Loff;}

        //! get first global j index (y-direction) on the subdomain
        //! (including ghost-nodes)
        inline int FirstJ(){return Moff;}

        //! get first global i index (x-direction) on the subdomain
        //! (including ghost-nodes)
        inline int FirstI(){return Noff;}

        //! get first global k index (z-direction) on the subdomain
        //! (excluding ghost-nodes)
        inline int FirstRealK(){return Loff0;}

        //! get first global j index (y-direction) on the subdomain
        //! (excluding ghost-nodes)
        inline int FirstRealJ(){return Moff0;}

        //! get first global i index (x-direction) on the subdomain
        //! (excluding ghost-nodes)
        inline int FirstRealI(){return Noff0;}

        //! get last global k index (z-direction) on the subdomain
        //! (including ghost-nodes)
        inline int LastK(){return Loff+lloc-1;}

        //! get last global j index (y-direction) on the subdomain
        //! (including ghost-nodes)
        inline int LastJ(){return Moff+mloc-1;}

        //! get last global i index (x-direction) on the subdomain
        //! (including ghost-nodes)
        inline int LastI(){return Noff+nloc-1;}

        //! get last global k index (z-direction) on the subdomain
        //! (excluding ghost-nodes)
        inline int LastRealK(){return Loff0+lloc0-1;}

        //! get last global j index (y-direction) on the subdomain
        //! (excluding ghost-nodes)
        inline int LastRealJ(){return Moff0+mloc0-1;}

        //! get last global i index (x-direction) on the subdomain
        //! (excluding ghost-nodes)
        inline int LastRealI(){return Noff0+nloc0-1;}

        //! get comm object
        inline Teuchos::RCP<Epetra_Comm> GetComm() { return comm;}

        //! Creates a communicator consisting of a 1D cut through the 3D comm.
        //! Will create a communicator containing all processes whose ranks
        //! differ only in dimension dim. For instance, GetProcRow(0)
        //! will return a comm with all subdomains with the same [ymin ymax]x[zmin zmax]
        //! and so al
        //! low communication in dimension 0 (x)
        Teuchos::RCP<Epetra_Comm> GetProcRow(int dim);

        //! The assembly map is created during the Decomp2D call.
        //! It contains ghost nodes between subdomains and should
        //! be used for interfacing with THCM (i.e. to compute
        //! global indices for the output of THCM on the subdomain)
        Teuchos::RCP<Epetra_Map> GetAssemblyMap(){return AssemblyMap;}

        //! The assembly surface map is a depth-averaged,
        //! single-unknown variant of the assembly map.
        Teuchos::RCP<Epetra_Map> GetAssemblySurfaceMap(){return AssemblySurfaceMap;}

        //! The standard map is created during the Decomp2D call.
        //! It does not contain ghost nodes and is therefore
        //! more suitable for the solve phase in Trilinos. The
        //! indexing scheme and decomposition are the same as in
        //! the assembly map.
        Teuchos::RCP<Epetra_Map> GetStandardMap(){return StandardMap;}

        //! The standard surface map is a depth-averaged,
        //! single-unknown variant of the standard map.
        Teuchos::RCP<Epetra_Map> GetStandardSurfaceMap(){return StandardSurfaceMap;}

        //! if you register a load-balancing object with this
        //! class, this map is made to optimize performance of
        //! linear solvers. The decomposition is still 2D, but
        //! the subdomains may no longer be rectangular.
        Teuchos::RCP<Epetra_Map> GetSolveMap(){return SolveMap;}

        //! The column map is created when the constructor is
        //! called. It is identical on all processes and contains
        //! all nodes on every subdomain. It is typically used as
        //! column map of a matrix, but might also serve for
        //! 'Allgather/Scatter' operations.
        Teuchos::RCP<Epetra_Map> GetColMap(){return ColMap;}

        //! Obtain local grid
        Grid GetLocalGrid() {return gridLoc_;}

        //! Obtain global grid
        const Grid& GetGlobalGrid() {return gridGlb_;}

        double Hdim() const {return gridGlb_.hdim_;}

        //@{ get physical bound of the subdomain
        double XminLoc() const {return gridLoc_.xmin_;}
        double XmaxLoc() const {return gridLoc_.xmax_;}
        double YminLoc() const {return gridLoc_.ymin_;}
        double YmaxLoc() const {return gridLoc_.ymax_;}
        //@}

        //@{ get physical bound of the global domain
        double Xmin() const {return gridGlb_.xmin_;}
        double Xmax() const {return gridGlb_.xmax_;}
        double Ymin() const {return gridGlb_.ymin_;}
        double Ymax() const {return gridGlb_.ymax_;}
        double Zmin() const {return gridGlb_.zmin_;}
        double Zmax() const {return gridGlb_.zmax_;}
        //@}

        //! returns true if the given local index 'li' is on a
        //! ghost-node. 'nun' is specified to indicate the unknowns
        //! per grid point
        bool IsGhost(int li, int nun_) const;
        bool IsGhost(int li) const { return IsGhost(li, dof_); };

        //! if true, the domain is periodic in the x-direction.
        bool IsPeriodic() const {return periodic;}

        //! create a 'solve' map with given number of unknowns per node.
        //! if depth_av=true, only one layer of grid points will
        //! be present. This function is not intended for creating submaps
        //! of the global map, see class SplitMatrix for that purpose.
        Teuchos::RCP<Epetra_Map> CreateStandardMap(int nun_, bool depth_av_=false) const;

        //! create a 'solve' map with given number of unknowns per node.
        //! if depth_av=true, only one layer of grid points will
        //! be present. This function is not intended for creating submaps
        //! of the global map, see class SplitMatrix for that purpose.
        Teuchos::RCP<Epetra_Map> CreateSolveMap(int nun_, bool depth_av_=false) const;

        //! create a 'assembly' map with given number of unknowns per node.
        //! if depth_av=true, only one layer of grid points will
        //! be present. This function is not intended for creating submaps
        //! of the global map, see class SplitMatrix for that purpose.
        Teuchos::RCP<Epetra_Map> CreateAssemblyMap(int nun_, bool depth_av_=false) const;

        //! this used to be a feature in trilinos_thcm, where it was possible
        //! to have an additional import operation between a 'standard' and a
        //! 'solve' map, so that subdomains with many land cells could be made
        //! bigger. Right now the standard and solve maps are the same, so this
        //! function always returns false.
        bool UseLoadBalancing() const {return false;}

        //@{ \name Data Transfer functions between the three map-types
        int Assembly2Standard(const Epetra_Vector& source, Epetra_Vector& target) const;

        int Assembly2StandardSurface(const Epetra_Vector& source,
                                     Epetra_Vector& target) const;

        int Assembly2Solve(const Epetra_Vector& source, Epetra_Vector& target) const;
        int Standard2Assembly(const Epetra_Vector& source, Epetra_Vector& target) const;

        int Standard2AssemblySurface(const Epetra_Vector& source,
                                     Epetra_Vector& target) const;

        int Standard2Solve(const Epetra_Vector& source, Epetra_Vector& target) const;
        int Solve2Assembly(const Epetra_Vector& source, Epetra_Vector& target) const;
        int Solve2Standard(const Epetra_Vector& source, Epetra_Vector& target) const;
        //@}
        //! we also offer this option for matrices, the others are not so important
        int Standard2Solve(const Epetra_CrsMatrix& source, Epetra_CrsMatrix& target) const;

    protected:

        //! communicator object
        Teuchos::RCP<Epetra_Comm> comm;

        //! see GetAssemblyMap() for a description
        Teuchos::RCP<Epetra_Map> AssemblyMap;

        //! see GetStandardMap() for a description
        Teuchos::RCP<Epetra_Map> StandardMap;

        //! see GetAssemblySurfaceMap() for a description
        Teuchos::RCP<Epetra_Map> AssemblySurfaceMap;

        //! see GetStandardSurfaceMap() for a description
        Teuchos::RCP<Epetra_Map> StandardSurfaceMap;

        //! see GetSolveMap() for a description --> deprecated
        Teuchos::RCP<Epetra_Map> SolveMap;

        //! see GetColMap() for a description
        Teuchos::RCP<Epetra_Map> ColMap;

        //! objects to transform the three vector types into one another
        Teuchos::RCP<Epetra_Import> as2std,std2sol,as2std_surf;

        int n,m,l; //!dimension of global domain
        int lloc,mloc,nloc; //! dimension of local subdomain (incl. ghost-nodes)
        int lloc0,mloc0,nloc0; //! dimension of local subdomain (excl. ghost-nodes)

        int Loff,Moff,Noff; //! offsets where the subdomain starts (with ghost nodes)
        int Loff0,Moff0,Noff0; //! offsets where the subdomain starts (without ghost nodes)

        //! position in processor array
        int pidL,pidM,pidN;

        //! size of processor array
        int npL,npM,npN;

        //! flag indicating if there is more than one subdomain in the x-direction
        bool xparallel;

        //! if true, the domain is periodic in the x-direction.
        //! This means that proc nx-1 has right neighbour proc
        //! 0 etc. and is useful for global computations.
        bool periodic;

        //! number of "ghost layers" (constant)
        /*! This is the number of cell layers imported
          from a neighbouring subdomains. The original
          number of grid points in each direction where
          there is a neighbour is increased by this amount
        */
        static const int numGhosts = 2;

        //! number of unknowns per grid cell
        int dof_;

        //! Number of auxiliary unknowns. Auxiliary unknowns do not
        //! have grid coordinates and are appended at the end of the
        //! ordinary map. Such rows can be used to compute
        //! precipitation.
        int aux_;

        //! Full grid representations
        Grid gridLoc_;
        const Grid gridGlb_;

    protected:

        void CommonSetup();

    private:

        //! private map generating function
        Teuchos::RCP<Epetra_Map> CreateMap(int noff_, int moff_, int loff_,
                                           int nloc_, int mloc_, int lloc_,
                                           int nun_,  bool addAux = false) const;

    };

}// namespace TRIOS

#endif // TRIOS_DOMAIN
