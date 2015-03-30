/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include "Teuchos_StandardCatchMacros.hpp"

#include "globdefs.H"

#include "TRIOS_Domain.H"

#include "OceanOutputXdmf.H"
#include "THCM.H"
#include "OceanGrid.H"
#ifdef HAVE_XDMF
#include <hdf5.h>
#endif



// write transformed velocity (otherwise you have to transform it
// yourself, only a transformation matrix is stored)
#define CARTESIAN_VELOCITY

using namespace Teuchos;
    
    //! constructor: opens the hdf5 file and stores grid if output==true.
    OceanOutputXdmf::OceanOutputXdmf(Teuchos::RCP<TRIOS::Domain> domain, 
           Teuchos::ParameterList& params, bool output) :
           domain_(domain)
      {
      counter_ = 0;
      comm_ = domain_->GetComm();
#ifdef HAVE_XDMF
    if (output)
      {
      hdf5_ = Teuchos::rcp(new EpetraExt::HDF5(*comm_));
      }
#endif
      xml_Domain_=Teuchos::rcp(new Teuchos::XMLObject("Domain"));
      xml_Grid1_=Teuchos::rcp(new Teuchos::XMLObject("Grid"));

      xdmf_version_=params.get("Xdmf Version","2.0");
      filename_=params.get("File Name","Trilinos_THCM");
      store_cartesian_=params.get("Store X/Y/Z Mesh",false);
      if (store_cartesian_)
        {
        radius_ = params.get("Earth-Radius for Visualization",10.0);
        xml_Grid2_=Teuchos::rcp(new Teuchos::XMLObject("Grid"));
        }
      write_uvw_=params.get("Store Velocity",true);
      write_p_=params.get("Store Pressure",true);
      write_T_=params.get("Store Temperature",true);
      write_S_=params.get("Store Salinity",true);
#ifdef HAVE_XDMF
    if (output)
      {
      hdf5_->Create((filename_+".h5").c_str());
      }
#endif
      
      grid_=Teuchos::rcp(new OceanGrid(domain_));
      cell_map_ = domain->CreateStandardMap(1,false);
      cell_vmap_ = domain->CreateStandardMap(3,false);

      scalarData_=Teuchos::rcp(new Epetra_Vector(*cell_map_));
      vectorData_=Teuchos::rcp(new Epetra_Vector(*cell_vmap_));
      
      // create the reference grid object
      Teuchos::RCP<Epetra_Vector> xc=grid_->GetXc(); 
      Teuchos::RCP<Epetra_Vector> yc=grid_->GetYc(); 
      Teuchos::RCP<Epetra_Vector> zc=grid_->GetZc(); 
      Teuchos::RCP<Epetra_Vector> xu=grid_->GetXu(); 
      Teuchos::RCP<Epetra_Vector> yv=grid_->GetYv(); 
      Teuchos::RCP<Epetra_Vector> zw=grid_->GetZw(); 
      if (output)
        {
        this->StoreGrid(xc,yc,zc,xu,yv,zw);
        }
        
      //this->Finalize(); 
      }
    
    //! destructor
    OceanOutputXdmf::~OceanOutputXdmf()
      {
      this->Finalize();
      }


  void OceanOutputXdmf::Finalize()
    {
#ifdef HAVE_XDMF
  if (hdf5_!=null)
    {
    INFO("Close HDF5 File...");
    hdf5_->Close();
    comm_->Barrier();
    }
#endif    
    }
    
  //! store the grid (has to be done only once)
  int OceanOutputXdmf::StoreGrid(Teuchos::RCP<Epetra_Vector> _xc, Teuchos::RCP<Epetra_Vector> _yc, Teuchos::RCP<Epetra_Vector> _zc, 
                                 Teuchos::RCP<Epetra_Vector> _xu,  Teuchos::RCP<Epetra_Vector> _yv,  Teuchos::RCP<Epetra_Vector> _zw)
    {
    Teuchos::RCP<Epetra_Vector> xc=_xc;
    Teuchos::RCP<Epetra_Vector> yc=_yc;
    Teuchos::RCP<Epetra_Vector> zc=_zc;
    Teuchos::RCP<Epetra_Vector> xu=_xu;
    Teuchos::RCP<Epetra_Vector> yv=_yv;
    Teuchos::RCP<Epetra_Vector> zw=_zw;
    
    
#ifdef HAVE_XDMF
    // write XML (light data)
    std::string domain_name = THCM::Instance().Label();
    xml_Domain_->addAttribute("Name",domain_name);

    int nglob = domain_->GlobalN();
    int mglob = domain_->GlobalM();
    int lglob = domain_->GlobalL();
    
    if (xdmf_version_=="1.0")
      {
      Error("Xdmf Version 1 output is not implemented!",__FILE__,__LINE__);
      }
    else
      {
      // create topology and geometry inside the reference grid
      // we put the nodes at the cell-centers of our mesh and make
      // all quantities node-centered, because cell-centered velocities
      // don't work in ParaView (I think)
      xml_Grid1_->addAttribute("GridType","Uniform");
      if (store_cartesian_)
        {
        xml_Grid2_->addAttribute("GridType","Uniform");
        }

      XMLObject topo1("Topology"); // spherical coordinates
      XMLObject topo2("Topology"); // cartesian coordinates

     std::string dims = toString(lglob) + " " + toString(mglob) + " " + toString(nglob);
        
      topo1.addAttribute("Dimensions", dims);
      topo2.addAttribute("Dimensions", dims);

      double hdim = THCM::Instance().hDim();
      double xmin = THCM::Instance().xMin();
      double xmax = THCM::Instance().xMax();
      double ymin = THCM::Instance().yMin();
      double ymax = THCM::Instance().yMax();
   
      if (store_cartesian_==true)
        {
        topo2.addAttribute("TopologyType","3DSMesh");
        xml_Grid2_->addChild(topo2);
        // cartesian coordinates
        double xpos,ypos,zpos,phi,lambda,z;
        double r=radius_;
        int pos=0;
        for (int k=0; k<zc->MyLength(); k++)
          for (int j=0; j<yc->MyLength(); j++)
            for (int i=0; i<xc->MyLength(); i++)
              {
              lambda = (*xc)[i];
              phi = (*yc)[j];
              z = (*zc)[k];
                
              xpos = (r+z)*cos(lambda)*cos(phi);
              ypos = (r+z)*sin(lambda)*cos(phi);
              zpos = (r+z)*sin(phi);

              (*vectorData_)[pos++]=xpos;
              (*vectorData_)[pos++]=ypos;
              (*vectorData_)[pos++]=zpos;
              }
        }
        
     
      // spherical coords are always included:
      topo1.addAttribute("TopologyType","3DRectMesh");
      xc = Teuchos::rcp(new Epetra_Vector(*_xc));
      yc = Teuchos::rcp(new Epetra_Vector(*_yc));
      zc = Teuchos::rcp(new Epetra_Vector(*_zc));
      xc->Scale(180.0/PI_);
      yc->Scale(180.0/PI_);
      //zc->Scale(-1.0);
        
      xml_Grid1_->addChild(topo1);
      XMLObject geom1("Geometry");
      XMLObject geom2("Geometry");
        XMLObject xyzdata("DataItem");
        XMLObject xdata("DataItem");
        XMLObject ydata("DataItem");
        XMLObject zdata("DataItem");
        if (store_cartesian_==true)
          {
          geom2.addAttribute("Type","XYZ");
          xyzdata.addAttribute("Name","XYZ");
          xyzdata.addAttribute("Format","HDF");
          xyzdata.addAttribute("NumberType","Float");
          xyzdata.addAttribute("Precision","8");
          xyzdata.addAttribute("Dimensions",toString(nglob*mglob*lglob)+" 3");
          xyzdata.addContent(filename_+".h5:/grid_xyz/Values");
          geom2.addChild(xyzdata);
          xml_Grid2_->addChild(geom2);
          }
        geom1.addAttribute("Type","VXVYVZ");
        xdata.addAttribute("Name","X");
        xdata.addAttribute("Format","HDF");
        xdata.addAttribute("NumberType","Float");
        xdata.addAttribute("Precision","8");
        xdata.addAttribute("Dimensions",toString(nglob));
        xdata.addContent(filename_+".h5:/grid_xc/Values");        
        geom1.addChild(xdata);                
        ydata.addAttribute("Name","Y");
        ydata.addAttribute("Format","HDF");
        ydata.addAttribute("NumberType","Float");
        ydata.addAttribute("Precision","8");
        ydata.addAttribute("Dimensions",toString(mglob));
        ydata.addContent(filename_+".h5:/grid_yc/Values");        
        geom1.addChild(ydata);                
        zdata.addAttribute("Name","Z");
        zdata.addAttribute("Format","HDF");
        zdata.addAttribute("NumberType","Float");
        zdata.addAttribute("Precision","8");
        zdata.addAttribute("Dimensions",toString(lglob));
        zdata.addContent(filename_+".h5:/grid_zc/Values");        
        geom1.addChild(zdata);                
      xml_Grid1_->addChild(geom1);
      }

      // create land mask array
      int i0 = domain_->FirstRealI()-domain_->FirstI()+1;
      int j0 = domain_->FirstRealJ()-domain_->FirstJ()+1;
      int k0 = domain_->FirstRealK()-domain_->FirstK()+1;
        
      int i1 = domain_->LastRealI()-domain_->FirstI()+1;
      int j1 = domain_->LastRealJ()-domain_->FirstJ()+1;
      int k1 = domain_->LastRealK()-domain_->FirstK()+1;
      int pos=0;
      for (int k=k0; k<=k1; k++)
        for (int j=j0; j<=j1; j++)
          for (int i=i0; i<=i1; i++)
            {
            (*scalarData_)[pos++] = (double)(grid_->landm(i,j,k));
            }

      // store an Xdmf file called "topo.xmf" containing the land mask
      XMLObject grid1 = xml_Grid1_->deepCopy(); // copy information from reference grid
      grid1.addAttribute("Name","Topography (Spherical Coordinates)");
      addAttribute(grid1, "Land Mask", "Scalar", 
                     filename_+".h5:/grid_landm/Values");
      

      if (store_cartesian_)
        {
        XMLObject grid2 = xml_Grid2_->deepCopy(); // copy information from reference grid
        grid2.addAttribute("Name","Topography (Cartesian Coordinates)");
        addAttribute(grid2, "Land Mask", "Scalar", 
                     filename_+".h5:/grid_landm/Values");                     
        this->WriteXmf(grid1,grid2,filename_+"_topo.xmf");
        }
      else
        {
        this->WriteXmf(grid1,filename_+"_topo.xmf");
        }

      // write HDF5 (heavy data)

      INFO("HDF5: Write Grid Data");

      hdf5_->Write("grid_xc",*xc);
      hdf5_->Write("grid_yc",*yc);
      hdf5_->Write("grid_zc",*zc);


      hdf5_->Write("grid_xu",*xu);      
      hdf5_->Write("grid_yv",*yv);
      hdf5_->Write("grid_zw",*zw);

      hdf5_->Write("grid_landm",*scalarData_);
      
      if (store_cartesian_)
        {
        hdf5_->Write("grid_xyz",(*vectorData_));      
        }      


      // the grid gets an attribute 'Vector Transformation'
      // which can be used in ParaView to transform the ve-
      // locity field to cartesian coordinates
      if (store_cartesian_)
        {
        //these vectors are such that u_cart = dot(transU,[u,v,w]) etc
        Epetra_Vector transU = *vectorData_;
        Epetra_Vector transV = *vectorData_;
        Epetra_Vector transW = *vectorData_;
        int pos=0;
        double lambda,phi,r;
        for (int k=k0; k<=k1; k++)
          {
          for (int j=j0; j<=j1; j++)
            {
            for (int i=i0; i<=i1; i++)
              {
              // to transform u/v/w to cartesian vectors
              lambda = (*xc)[i-1]; // adjust to 0-based local index
              phi = (*yc)[j-1];
              r = radius_+(*zc)[k-1];

              transU[pos]   =-r*sin(lambda)*cos(phi);
              transU[pos+1] =-r*cos(lambda)*sin(phi);
              transU[pos+2] = cos(lambda)*cos(phi);
              transV[pos]   = r*cos(lambda)*cos(phi);
              transV[pos+1] =-r*sin(lambda)*sin(phi);
              transV[pos+2] = sin(lambda)*cos(phi);
              transW[pos]   = 0.0;
              transW[pos+1] = r*cos(phi);
              transW[pos+2] = sin(phi);
              pos+=3;
              }
            }
          }
        hdf5_->Write("sph2cartu",transU);
        hdf5_->Write("sph2cartv",transV);
        hdf5_->Write("sph2cartw",transW);
        }

/////////////////////////////////

#endif
      return 0;      
      }


    
  //! add solution at new time-step to file
  int OceanOutputXdmf::Store(const Epetra_Vector& sol, double t,bool xmf_out)
    {
    // (a) import to grid
    grid_->ImportData(sol);
    if (xmf_out) 
      {
#ifdef HAVE_XDMF
      counter_++;
      std::stringstream ss;
      ss<<std::setw(4)<<std::setfill('0')<<counter_;
     std::string number = ss.str();
     std::string groupname="state"+number;
      
      INFO("Saving Xdmf "<<groupname<<" at t=" << t);

    XMLObject grid = xml_Grid1_->deepCopy(); // copy information from reference grid
    FillXmlGrid(grid,"Spherical Coordinates",groupname,t);
    WriteXmf(grid, "thcm_state_"+number+".xmf");
    
    if (store_cartesian_)
      {
      XMLObject grid2 = xml_Grid2_->deepCopy(); // copy information from reference grid
      FillXmlGrid(grid2,"Cartesian Coordinates",groupname,t,true);
      WriteXmf(grid,grid2, "thcm_state_"+number+".xmf");
      }

    // get the solution components as parallel Epetra Vectors:
    
    //write heavy data

//      hdf5_->Write("/",groupname+"_index",counter_);
      
//      hdf5_->Write("/",groupname+"_time",t);

      // relevant index range: note that the OceanGrid class contains ghost-nodes
      // between subdomains which we do not want to include in the output.
      // the range of interest is from 1 to n whereas THCM arrays go from 0 to n+1
      // (except for the u and v arrays which go from 0 to n)
      
      // note: these ive 1-based indexing as used by the grid object
      int i0 = domain_->FirstRealI()-domain_->FirstI()+1;
      int j0 = domain_->FirstRealJ()-domain_->FirstJ()+1;
      int k0 = domain_->FirstRealK()-domain_->FirstK()+1;
        
      int i1 = domain_->LastRealI()-domain_->FirstI()+1;
      int j1 = domain_->LastRealJ()-domain_->FirstJ()+1;
      int k1 = domain_->LastRealK()-domain_->FirstK()+1;
        
      if (write_uvw_)
        {
        // construct a vector with u/v/w interpolated to the cell-centers
        // and in interleaved ordering:
        int pos=0;
        double lambda,phi,r;
        double utmp,vtmp,wtmp;
        Teuchos::RCP<Epetra_Vector> xc,yc,zc;
        if (store_cartesian_)
          {
          xc=grid_->GetXc();
          yc=grid_->GetYc();
          zc=grid_->GetZc();
          }
        // compute velocity at cell centers (interpolation)
        for (int k=k0; k<=k1; k++)
          {
          for (int j=j0; j<=j1; j++)
            {
            for (int i=i0; i<=i1; i++)
              {
              (*vectorData_)[pos++]=0.25*(grid_->u(i,j,k)   + grid_->u(i-1,j,k)
                                    + grid_->u(i,j-1,k) + grid_->u(i-1,j-1,k));
              (*vectorData_)[pos++]=0.25*(grid_->v(i,j,k)   + grid_->v(i-1,j,k)
                                    + grid_->v(i,j-1,k) + grid_->v(i-1,j-1,k));
              (*vectorData_)[pos++]=0.5*(grid_->w(i,j,k)   + grid_->w(i,j,k-1));
#ifdef CARTESIAN_VELOCITY
              if (store_cartesian_)
                {
                // transform u/v/w to cartesian vectors
                lambda = (*xc)[i-1]; // adjust to 0-based local index
                phi = (*yc)[j-1];
                r = radius_+(*zc)[k-1];
                // multiply velocity by jacobian d(x,y,z)/d(lambda,phi,z')
                utmp = -r*sin(lambda)*cos(phi)*(*vectorData_)[pos-3]
                      -r*cos(lambda)*sin(phi)*(*vectorData_)[pos-2]
                      +  cos(lambda)*cos(phi)*(*vectorData_)[pos-1];
                vtmp =  r*cos(lambda)*cos(phi)*(*vectorData_)[pos-3]
                      -r*sin(lambda)*sin(phi)*(*vectorData_)[pos-2]
                      +  sin(lambda)*cos(phi)*(*vectorData_)[pos-1];
                wtmp =  r*cos(phi)*(*vectorData_)[pos-2]
                      +  sin(phi)*(*vectorData_)[pos-1];
                (*vectorData_)[pos-3]=utmp;
                (*vectorData_)[pos-2]=vtmp;
                (*vectorData_)[pos-1]=wtmp;
                }
#endif          
              }
            }            
          }
        hdf5_->Write(groupname+"_velocity",(*vectorData_));
        }
      if (write_p_)
        {
        int pos=0;
        for (int k=k0; k<=k1; k++)
          for (int j=j0; j<=j1; j++)
            for (int i=i0; i<=i1; i++)
              {
              (*scalarData_)[pos++]=grid_->p(i,j,k);
              }       
        hdf5_->Write(groupname+"_pressure",(*scalarData_)); 
        }
      if (write_T_)
        {
        int pos=0;
        for (int k=k0; k<=k1; k++)
          for (int j=j0; j<=j1; j++)
            for (int i=i0; i<=i1; i++)
              {
              (*scalarData_)[pos++]=grid_->T(i,j,k);
              }       
        hdf5_->Write(groupname+"_temperature",(*scalarData_));        
        }
      if (write_S_)
        {
        int pos=0;
        for (int k=k0; k<=k1; k++)
          for (int j=j0; j<=j1; j++)
            for (int i=i0; i<=i1; i++)
              {
              (*scalarData_)[pos++]=grid_->S(i,j,k);
              }       
        hdf5_->Write(groupname+"_salinity",(*scalarData_));        
        }
#else
  INFO("Warning: Xdmf output is disabled, define HAVE_XDMF to change this");
  INFO ("("<<__FILE__<<", line "<<__LINE__<<")");
#endif
      }
      return 0;
      }

void OceanOutputXdmf::addAttribute(XMLObject& state,std::string name, 
                            std::string type,std::string hdf5location)
  {
  int N = domain_->GlobalN();
  int M = domain_->GlobalM();
  int L = domain_->GlobalL();
 std::string dimstring = toString(L)+" "+toString(M)+" "+toString(N);
  XMLObject attr("Attribute");
    attr.addAttribute("Name",name);
    // note: we put everything at nodes which we put at our cell-centers.
    //       So we need to interpolate u/v/w to the center
    attr.addAttribute("Center","Node");
    attr.addAttribute("Type",type);
    if (xdmf_version_=="1.0")
      {
      XMLObject data("DataStructure");
      data.addAttribute("Format","HDF");
      data.addAttribute("DataType","Float");
      data.addAttribute("Precision","8");
//    data.addAttribute("Dimensions",toString(imaxus_)+" "+toString(jmaxus_));
      if (type=="Scalar")
        {
        data.addAttribute("Dimensions",dimstring);
        }
      else if (type=="Vector")
       {
       data.addAttribute("Dimensions",dimstring+" 3");
       }
     else
       {
       Error("Invalid data type: "+type,__FILE__,__LINE__);
       }
      data.addContent(hdf5location);
      attr.addChild(data);
      }
    else
      {
      XMLObject data("DataItem");
      data.addAttribute("ItemType","Uniform");
      data.addAttribute("Format","HDF");
      data.addAttribute("NumberType","Float");
      data.addAttribute("Precision","8");
      if (type=="Scalar")
        {
        data.addAttribute("Dimensions",dimstring);
        }
      else if (type=="Vector")
       {
       data.addAttribute("Dimensions",dimstring+" 3");
       }
     else
       {
       Error("Invalid data type: "+type,__FILE__,__LINE__);
       }
      data.addContent(hdf5location);
      attr.addChild(data);
      }
  state.addChild(attr);
  }

void OceanOutputXdmf::FillXmlGrid(XMLObject& grid,
           std::string gridname,std::string groupname,
            double t, bool uvw_transform)
  {
  grid.addAttribute("Name",gridname);
  XMLObject time("Time");
  time.addAttribute("Value",toString(t));
  grid.addChild(time);
  //write light-weight data
  if (write_uvw_)
    {
    addAttribute(grid, "Velocity", "Vector", 
                 filename_+".h5:/"+groupname+"_velocity/Values");
    if (uvw_transform)
      {
      addAttribute(grid, "Sph2CartU", "Vector",
                 filename_+".h5:/"+groupname+"_sph2cartu/Values");
      addAttribute(grid, "Sph2CartV", "Vector",
                 filename_+".h5:/"+groupname+"_sph2cartv/Values");
      addAttribute(grid, "Sph2CartW", "Vector",
                 filename_+".h5:/"+groupname+"_sph2cartw/Values");
      }
    }
  if (write_p_)
    {
    addAttribute(grid, "Pressure", "Scalar",
                filename_+".h5:/"+groupname+"_pressure/Values");
    }
  if (write_T_)
    {
    addAttribute(grid, "Temperature", "Scalar",
                 filename_+".h5:/"+groupname+"_temperature/Values");
    }
  if (write_S_)
    {
    addAttribute(grid, "Salinity", "Scalar",
                 filename_+".h5:/"+groupname+"_salinity/Values");
    }
  }

// store single-grid XML file
void OceanOutputXdmf::WriteXmf(XMLObject& grid,std::string filename)
  {
  if (comm_->MyPID()==0)
    {
    std::ofstream xmlfile(filename.c_str());
    xmlfile<<"<?xml version=\"1.0\" ?>\n";
    xmlfile<<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";

    XMLObject xdmf("Xdmf");
    xdmf.addAttribute("Version",xdmf_version_);
    
    XMLObject domain = xml_Domain_->deepCopy();
    domain.addChild(grid);
    xdmf.addChild(domain);
    xmlfile<<xdmf;
    xmlfile.close();
    }
  }

// store two-grid XML file
void OceanOutputXdmf::WriteXmf(XMLObject& grid1,XMLObject& grid2,std::string filename)
  {
  if (comm_->MyPID()==0)
    {
    std::ofstream xmlfile(filename.c_str());
    xmlfile<<"<?xml version=\"1.0\" ?>\n";
    xmlfile<<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";

    XMLObject xdmf("Xdmf");
    xdmf.addAttribute("Version",xdmf_version_);
    
    XMLObject domain = xml_Domain_->deepCopy();
    domain.addChild(grid1);
    domain.addChild(grid2);
    xdmf.addChild(domain);
    xmlfile<<xdmf;
    xmlfile.close();
    }
  }
