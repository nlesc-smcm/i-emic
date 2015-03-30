/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#ifndef OCEAN_OUTPUT_XDMF_H
#define OCEAN_OUTPUT_XDMF_H


#include "globdefs.H"

// in this file we make an exception and use all of Teuchos and std,
// because EpetraExt_HDF5.hpp can otherwise not be included.
using namespace std;
using namespace Teuchos;

#ifdef HAVE_XDMF
#include "EpetraExt_HDF5.h"
#endif
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Array.hpp"
#include "Epetra_Comm.h"


namespace Teuchos {
  class ParameterList;
}

class Epetra_Vector;
namespace TRIOS {
class Domain;
}
class OceanGrid;

//!                                                              
//! an interface to the Trilinos HDF5 and XML implementations    
//! to allow storing series of solutions in the Xdmf format.     
//! Xdmf (eXtensible Data Model and Format, www.xdmf.org) stores 
//! 'light-weight data' in an XML file (ascii) and 'heavy data'  
//! in an HDF5 (Hierarchical Data Format, binary) file. The Xdmf 
//! format is supported by popular postprocessing and visualiza- 
//! tion tools such as ParaView, IDL and Ensight.                
//! The HDF5 file can be read and written in parallel, which     
//! makes it the format of choice in Trilinos-THCM.              
//!                                                              
class OceanOutputXdmf
  {
  public:
    
    //! constructor: if output==true, the HDF5 file is opened and the grid/topography is 
    //! stored. Otherwise, only an OceanGrid object is created.
    OceanOutputXdmf(Teuchos::RCP<TRIOS::Domain> domain, 
                    Teuchos::ParameterList& params, bool output=true);
    
    //! destructor
    ~OceanOutputXdmf();
    
    //! add current solution to file
    int Store(const Epetra_Vector& solution, double t,bool xmf_out);
    
    //! return grid object
    RCP<OceanGrid> getGrid() { return grid_; }

  protected:
        
    //! filename
    std::string filename_;
    
    //! grid decomposition object
    RCP<TRIOS::Domain> domain_;
    
    //! grid representation of solution vectors
    RCP<OceanGrid> grid_;

    //! communicator
    RCP<Epetra_Comm> comm_;
    
    //! map for single cell-centered quantity (p/T/S)
    RCP<Epetra_Map> cell_map_;

    //! map for vector cell-centered quantity 
    //! (interleaved velocity [uvw] mapped to cell centers)
    RCP<Epetra_Map> cell_vmap_;
    
    //! which components are to be written?
    bool write_uvw_,write_p_,write_T_,write_S_;

#ifdef HAVE_XDMF
    //! This is the HDF5 file manager for the heavy data
    RCP<EpetraExt::HDF5> hdf5_;
#endif
    
    //! This is the XML object for the light-weight data
    RCP<XMLObject> xml_Grid1_, xml_Grid2_, xml_Domain_;
    
    //! number of plots so far
    int counter_;
    
    //! which version of Xdmf do we write?
    std::string xdmf_version_;
    
    //! store 'real world' coordinates in addition to polar ones 
    bool store_cartesian_;
    
    //! if we store x/y/z coordinates, this is the radius of the planet for visualization:
    double radius_;
    
    //! vectors for parallel output
    RCP<Epetra_Vector> scalarData_,vectorData_;

  private:
  
    //! write light-weight data (xml file) at end of simpulation and close hdf5 file
    //! This function is called automatically by the destructor
    void Finalize();
  
    //! store the grid (has to be done only once)
    //! the xi/yj/zk arguments determine the cell centers, whereas the
    //! x/y/z arguments determine the cell corners (nodes)
    int StoreGrid(RCP<Epetra_Vector> xi, RCP<Epetra_Vector> yj, RCP<Epetra_Vector> zk, RCP<Epetra_Vector> x, RCP<Epetra_Vector> y, RCP<Epetra_Vector> z);
 
 
    //! private function for adding a variable to a grid
    void addAttribute(XMLObject& state, std::string name, 
                      std::string type, std::string hdf5location);   

    //! internal function to put attributes in XML grid entry  
    //! (set uvw_transform=true for cartesian coordinate grids)
    void FillXmlGrid(XMLObject& grid, std::string gridname, std::string groupname,
                     double t, bool uvw_transform=false);
    
    //! stores the liht-weight data in a single-grid file.
    void WriteXmf(XMLObject& grid, std::string filename);

    //! stores the liht-weight data in a two-grid file.
    void WriteXmf(XMLObject& grid1, XMLObject& grid2, std::string filename);
    
    
  };

#endif
