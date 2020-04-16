/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
/**********************************************************************
 * Modified and extended by Erik, Utrecht University 2014/15/16/17    *
 * contact: t.e.mulder@uu.nl                                          *
 **********************************************************************/

#include "Utils.H"

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_HDF5.h"

#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

#include "TRIOS_Domain.H"

#include "Combined_MultiVec.H"
#include "ComplexVector.H"

#include <functional> // for std::hash

using ConstIterator = Teuchos::ParameterList::ConstIterator;
//========================================================================================

//! simple ddot wrapper
double Utils::dot(std::vector<double> &vec1, std::vector<double> &vec2)
{
    assert(vec1.size() == vec2.size());
    int dim = (int) vec1.size();
    int incX = 1;
    int incY = 1;
    double dot;
    dot = ddot_(&dim, &vec1[0], &incX, &vec2[0], &incY);
    return dot;
}

//! simple summation
double Utils::sum(std::vector<double> &vec)
{
    double result = 0.0;
    for (auto &e: vec)
        result += e;
    return result;
}

//! count number of nonzeros using a threshold
int Utils::nnz(Teuchos::RCP<Epetra_Vector> vec, double t)
{
    int numMyElements = vec->Map().NumMyElements();
    int numInts = 0;
    int numIntsGlob = 0;
    for (int i = 0; i != numMyElements; ++i)
    {
        if (std::abs((*vec)[i]) > t)
            numInts++;
    }
    vec->Map().Comm().SumAll(&numInts, &numIntsGlob, 1);
    return numIntsGlob;
}

//! Clone an Epetra_Vector
Teuchos::RCP<Epetra_Vector> Utils::clone(
    Teuchos::RCP<const Epetra_Vector> const &vec)
{
    return Teuchos::rcp(new Epetra_Vector(*vec));
}

//! Clone a Combined_MultiVec
std::shared_ptr<Combined_MultiVec> Utils::clone(
    std::shared_ptr<const Combined_MultiVec> const &vec)
{
    std::vector<Teuchos::RCP<Epetra_Vector>> vecs;
    for (int i = 0; i < vec->Size(); i++)
        vecs.push_back(
            Teuchos::rcp(
                new Epetra_Vector(
                    *Teuchos::rcp_dynamic_cast<const Epetra_Vector>((*vec)(i)))));

    if (vec->Size() == 2)
        return std::make_shared<Combined_MultiVec>(vecs[0], vecs[1]);
    else if (vec->Size() == 3)
        return std::make_shared<Combined_MultiVec>(vecs[0], vecs[1], vecs[2]);

    return std::make_shared<Combined_MultiVec>(*vec);
}

//! Obtain 2-norm of std::vector<double> using ddot
double Utils::norm(std::vector<double> &vec)
{
    int dim = (int) vec.size();
    int incX = 1;
    int incY = 1;
    double dot;
    dot = ddot_(&dim, &vec[0], &incX, &vec[0], &incY);
    return sqrt(dot);
}

//! Obtain 2-norm of first vec in multivec, more convenient interface
double Utils::norm(Epetra_MultiVector &vec)
{
    std::vector<double> norm(vec.NumVectors());
    CHECK_ZERO(vec.Norm2(&norm[0]));
    return norm[0];
}

//! Obtain 2-norm of first vec in multivec, more convenient interface
double Utils::norm(Combined_MultiVec &vec)
{
    std::vector<double> norm(vec.NumVectors());
    CHECK_ZERO(vec.Norm2(&norm[0]));
    return norm[0];
}

//! Obtain first inf-norm of Epetra_MultiVector
double Utils::normInf(Epetra_MultiVector &vec)
{
    std::vector<double> norm(vec.NumVectors());
    CHECK_ZERO(vec.NormInf(&norm[0]));
    return norm[0];
}

//! Obtain first inf-norm of Combined_MultiVec
double Utils::normInf(Combined_MultiVec &vec)
{
    std::vector<double> norm(vec.NumVectors());
    CHECK_ZERO(vec.NormInf(&norm[0]));
    return norm[0];
}

//! Return the norm for a specific field using domain decomposition
//! given by <domain>. Field is specified by <unknown> (1-6 for an
//! ocean vector)
double Utils::normOfField(Teuchos::RCP<Epetra_MultiVector> vec,
                          Teuchos::RCP<TRIOS::Domain> domain, int unknown)
{
    Teuchos::RCP<Epetra_Map> rowMap  = domain->GetSolveMap();

    Teuchos::RCP<Epetra_Map> mapid   =
        Utils::CreateSubMap(*rowMap, domain->Dof(), unknown);

    Teuchos::RCP<Epetra_Import> importid =
        Teuchos::rcp(new Epetra_Import(*rowMap, *mapid));

    Epetra_Vector xid(*mapid);
    CHECK_ZERO(xid.Export(*vec, *importid, Zero));
    double norm;
    xid.Norm2(&norm);
    return norm;
}


//! Compute column sums. This is imitated from
//! Epetra_CrsMatrix::NormOne()
void Utils::colSums(Epetra_CrsMatrix const &mat, Epetra_Vector &sums)
{
    if (!mat.Filled())
    {
        ERROR("Matrix not filled", __FILE__, __LINE__);
    }

    if (!sums.Map().SameAs(mat.DomainMap()))
    {
        CHECK_ZERO(sums.ReplaceMap(mat.DomainMap()));
    }

    double *sumVals = (double *) sums.Values();
    Epetra_MultiVector *tmp = 0;
    int NumCols = mat.NumMyCols();


    if (mat.Importer() != 0) // non-trivial importer
    {
        tmp     = new Epetra_Vector(mat.ColMap()); // import vector
        sumVals = (double *) tmp->Values(); // point sumVals to tmp values
    }

    for (int i = 0; i < NumCols; ++i)
        sumVals[i] = 0.0;

    int     NumEntries;
    int    *ColIndices;
    double *RowValues;
    for (int i = 0; i < mat.NumMyRows(); ++i)
    {
        NumEntries = mat.NumMyEntries(i);
        CHECK_ZERO(mat.ExtractMyRowView(i, NumEntries, RowValues, ColIndices));
        for (int j = 0; j < NumEntries; ++j)
            sumVals[ColIndices[j]] += RowValues[j];
    }

    if (mat.Importer() != 0) // non-trivial importer
    {
        sums.PutScalar(0.0);
        CHECK_ZERO(sums.Export(*tmp, *mat.Importer(), Add));
    }

    if (tmp != 0)
        delete tmp;
}

//! Update std::vector<double>, result is stored in B
//! B = scalarA*A+scalarB*B
void Utils::update(double scalarA, std::vector<double> &A,
                   double scalarB, std::vector<double> &B)
{
    int N = (int) B.size();
    int incX = 1; int incY = 1;

    // scale B
    dscal_(&N, &scalarB, &B[0], &incX);

    // update B
    daxpy_(&N, &scalarA, &A[0], &incX, &B[0], &incY);
}

//! Print std vector
void Utils::print(std::vector<double> const &vec, std::string const &fname)
{
    std::ofstream file;
    file.open(fname);
    for (auto &el: vec)
        file << el << std::endl;
}

//! Print surface mask to INFO and file
void Utils::printSurfaceMask(Utils::MaskStruct const &mask)
{
    printSurfaceMask(mask.global_surface, mask.label, mask.n);
}

//! Print surface mask to INFO and file fname
void Utils::printSurfaceMask(std::shared_ptr<std::vector<int> > mask,
                             std::string const &fname,
                             int nDim)
{
    INFO("\nPrinting surface mask to " << fname);
    
    std::ostringstream string;
    std::ofstream smask;

    std::vector<std::string> stringvec;
    int ctr  = 0;
    int ctr0 = 0;
    int ctr1 = 0;
    for (auto &el: *mask)
    {
        ctr++;
        string << el;
        if (ctr % nDim == 0)
        {
            stringvec.push_back(string.str());
            string.str("");
            string.clear();
        }
        if (el == 0)
            ctr0++;
        else
            ctr1++;
    }
    INFO("   surface mask zeros: " << ctr0 << ", ones: " << ctr1 << '\n');
    // Reverse write to info and output file
    smask.open(fname);
    for (auto i = stringvec.rbegin(); i != stringvec.rend(); ++i)
    {
        INFO(i->c_str());
        smask << i->c_str() << '\n';
    }
    smask.close();
    INFO("");
}

//! We do a recursion for originalPars, not for dominantPars
void Utils::overwriteParameters(Teuchos::RCP<Teuchos::ParameterList> originalPars,
                                Teuchos::RCP<Teuchos::ParameterList> dominantPars)
{
    if (originalPars->begin() == originalPars->end() ||
        dominantPars->begin() == dominantPars->end() )
        WARNING("Feeding empty list into overwriteParameters", __FILE__, __LINE__);

    for (ConstIterator i = originalPars->begin(); i != originalPars->end(); ++i)
    {
        const Teuchos::ParameterEntry &entry_i = originalPars->entry(i);
        const std::string &name_i = originalPars->name(i);

        if (entry_i.isList()) // skipping the sublists first
        {
            Teuchos::RCP<Teuchos::ParameterList> sublist =
                Teuchos::rcp(&originalPars->sublist(name_i), false);

            std::stringstream ss;
            ss << originalPars->name() << "::" << name_i;
            sublist->setName(ss.str());
            overwriteParameters(sublist, dominantPars);
        }
        else
        {
            if (dominantPars->isParameter(name_i))
            {
                const Teuchos::ParameterEntry &entry_j =
                    dominantPars->getEntry(name_i);

                if (entry_i != entry_j)
                {

                    INFO(" " << originalPars->name() << "::"
                         << name_i << " = " << entry_i << " <-- "
                         << dominantPars->name() << "::" << name_i
                         << " = " << entry_j);
                    originalPars->setEntry(name_i, entry_j);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
Teuchos::RCP<Teuchos::ParameterList> Utils::obtainParams(std::string const &str,
                                                         std::string const &name)
{
    Teuchos::RCP<Teuchos::ParameterList> pars = rcp(new Teuchos::ParameterList);

    std::ifstream file(str);
    if (!file)
    {
        WARNING(str << ", " << name <<
                " not found, continuing with defaults at your own risk!",
                __FILE__, __LINE__);
    }
    else if (str == "dummy")
    {
        INFO("Continuing with dummy " << name);
    }
    else
    {
        Teuchos::updateParametersFromXmlFile(str.c_str(), pars.ptr());
    }
    pars->setName(name.c_str());
    return pars;
}

//-----------------------------------------------------------------------------
void Utils::obtainParams(Teuchos::RCP<Teuchos::ParameterList> pars,
                         std::string const &str,
                         std::string const &name)
{
    std::ifstream file(str);
    if (!file)
    {
        WARNING(str << ", " << name <<
                " not found, continuing with defaults at your own risk!",
                __FILE__, __LINE__);
    }
    else if (str == "dummy")
    {
        INFO("Continuing with dummy " << name);
    }
    else
    {
        Teuchos::updateParametersFromXmlFile(str.c_str(),
                                             Teuchos::sublist(pars, name).ptr());
    }
}

//=============================================================================
size_t Utils::hash(Teuchos::RCP<Epetra_MultiVector> vec)
{
    // Using an XOR and rotate hash on std::hash
    // This seems to work OK
    size_t seed = 2;
    std::hash<double> double_hash;

    // Loop over the vectors in the multivector, create hash
    int numVectors = vec->NumVectors();
    int numMyElements;
    for (int j = 0; j < numVectors; ++j)
    {
        numMyElements = (*vec)(j)->Map().NumMyElements();
        for (int i = 0; i < numMyElements; ++i)
            seed ^= double_hash((*(*vec)(0))[i]) + (seed << 6) + (seed >> 2);
    }

    return seed;
}


//============================================================================
void Utils::save(Teuchos::RCP<const Epetra_MultiVector> vec, std::string const &filename)
{
    INFO("Saving " << vec->Label() << " to " << filename);
    std::ostringstream fname;
    fname << filename << ".h5";
    EpetraExt::HDF5 HDF5(vec->Map().Comm());
    HDF5.Create(fname.str());

    // plot scripts will expect an entry called "State"
    HDF5.Write("State", *vec);
}

//============================================================================
void Utils::load(Teuchos::RCP<Epetra_MultiVector> vec, std::string const &filename)
{
    std::ostringstream fname;
    fname << filename;

    // Add extension if necessary
    std::string ext = ".h5";
    if (filename.length() < 3 ||
        !std::equal(ext.rbegin(), ext.rend(), filename.rbegin()))
        fname << ext;

    // Check whether file exists
    std::ifstream file(fname.str());
    if (!file)
    {
        WARNING("Can't open " << fname.str()
                << ", continue with trivial " << vec->Label(),
                __FILE__, __LINE__);

        // create trivial vector
        vec->PutScalar(0.0);
        return;
    }
    else file.close();

    INFO("Loading from " << fname.str() << " into " << vec->Label());

    EpetraExt::HDF5 HDF5(vec->Map().Comm());
    HDF5.Open(fname.str());
    Epetra_MultiVector *readState;

    // Check contents
    if (!HDF5.IsContained("State"))
    {
        ERROR("The group <State> is not contained in hdf5 " << filename,
              __FILE__, __LINE__);
    }

    HDF5.Read("State", readState);

    if ( readState->GlobalLength() != vec->Map().NumGlobalElements() )
    {
        ERROR("Incompatible number of elements", __FILE__, __LINE__);
    }

    // Create importer
    // target map: vector map
    // source map: state with linear map as read by HDF5.Read
    Teuchos::RCP<Epetra_Import> lin2solve =
        Teuchos::rcp(new Epetra_Import(vec->Map(),
                                       readState->Map()));

    // Import state from HDF5 into state_ datamember
    CHECK_ZERO(vec->Import(*((*readState)(0)), *lin2solve, Insert));

    delete readState;
}

//============================================================================
void Utils::save(std::shared_ptr<const Combined_MultiVec> vec, std::string const &filename)
{
    for (int i = 0; i != vec->Size(); ++i)
    {
        std::stringstream fname;
        fname << filename << "." << i;
        save( (*vec)(i), fname.str());
    }
}

//============================================================================
void Utils::load(std::shared_ptr<Combined_MultiVec> vec, std::string const &filename)
{
    for (int i = 0; i != vec->Size(); ++i)
    {
        std::stringstream fname;
        fname << filename << "." << i;
        load( (*vec)(i), fname.str());
    }
}

//============================================================================
void Utils::save(Combined_MultiVec const &vec, std::string const &filename)
{
    for (int i = 0; i != vec.Size(); ++i)
    {
        std::stringstream fname;
        fname << filename << "." << i;
        save( vec(i), fname.str());
    }
}

//============================================================================
// save eigenvectors based on combined_multivec
void Utils::saveEigenvectors(std::vector<ComplexVector<Combined_MultiVec> > const &eigvs,
                             std::vector<std::complex<double> > const &alpha,
                             std::vector<std::complex<double> > const &beta,
                             std::string const &filename)
{
    // Iterate over number of combined multivectors. Each multivector
    // gets its own HDF5 export process and corresponding file.
    int size = eigvs[0].real.Size();
    for (int i = 0; i != size; ++i)
    {
        std::stringstream ss, groupNameRe, groupNameIm;
        ss << filename << "." << i << ".h5";

        // Create HDF5 destination. We assume that the real and
        // imaginary part have the same map.
        EpetraExt::HDF5 HDF5(eigvs[0].real(i)->Map().Comm());

        HDF5.Create(ss.str().c_str());

        size_t ctr = 0;
        for (auto &vec: eigvs)
        {
            groupNameRe << "EV_Real_" << ctr;
            groupNameIm << "EV_Imag_" << ctr;

            HDF5.Write( groupNameRe.str().c_str(),
                        *vec.real(i) );
            HDF5.Write( groupNameIm.str().c_str(),
                        *vec.imag(i) );

            // clear stringstreams
            groupNameRe.str("");
            groupNameRe.clear();
            groupNameIm.str("");
            groupNameIm.clear();

            INFO("a / b = " << alpha[ctr] / beta[ctr]);
            ctr++;
        }

        // save the eigenvalues to all hdf5 files.
        saveEigenvalues(HDF5, alpha, beta, (int) eigvs.size());
    }
}

//============================================================================
// save eigenvectors based on epetra_vector
void Utils::saveEigenvectors(std::vector<ComplexVector<Epetra_Vector> > const &eigvs,
                             std::vector<std::complex<double> > const &alpha,
                             std::vector<std::complex<double> > const &beta,
                             std::string const &filename)
{
    std::stringstream ss, groupNameRe, groupNameIm;
    ss << filename << ".h5";

    // We assume the imaginary and real part of the ComplexVector have the same Map
    EpetraExt::HDF5 HDF5(eigvs[0].real.Map().Comm());

    HDF5.Create(ss.str().c_str());

    size_t ctr = 0;
    for (auto &vec: eigvs)
    {
        groupNameRe << "EV_Real_" << ctr;
        groupNameIm << "EV_Imag_" << ctr;

        HDF5.Write(groupNameRe.str().c_str(),
                   vec.real );
        HDF5.Write(groupNameIm.str().c_str(),
                   vec.imag );

        groupNameRe.str("");
        groupNameRe.clear();
        groupNameIm.str("");
        groupNameIm.clear();

        INFO("a / b = " << alpha[ctr] / beta[ctr]);
        ctr++;
    }

    saveEigenvalues(HDF5, alpha, beta, (int) eigvs.size());
}

//=============================================================================
void Utils::saveEigenvalues(EpetraExt::HDF5 &HDF5,
                            std::vector<std::complex<double> > const &alpha,
                            std::vector<std::complex<double> > const &beta,
                            int numEigs)
{
    std::stringstream nameRe, nameIm;
    HDF5.Write("MetaData", "NumEigs", numEigs);

    // Separate real and imaginary parts
    std::vector<double> alphaRe(numEigs, 0.0);
    std::vector<double> alphaIm(numEigs, 0.0);
    std::vector<double> betaRe(numEigs, 0.0);
    std::vector<double> betaIm(numEigs, 0.0);

    for (int i = 0; i != numEigs; ++i)
    {
        alphaRe[i] = alpha[i].real();
        alphaIm[i] = alpha[i].imag();
        betaRe[i]  = beta[i].real();
        betaIm[i]  = beta[i].imag();
    }

    HDF5.Write("EigenValues", "AlphaRe", H5T_NATIVE_DOUBLE, numEigs, &alphaRe[0]);
    HDF5.Write("EigenValues", "AlphaIm", H5T_NATIVE_DOUBLE, numEigs, &alphaIm[0]);
    HDF5.Write("EigenValues", "BetaRe",  H5T_NATIVE_DOUBLE, numEigs, &betaRe[0]);
    HDF5.Write("EigenValues", "BetaIm",  H5T_NATIVE_DOUBLE, numEigs, &betaIm[0]);
}

//=============================================================================
Teuchos::RCP<Epetra_Vector> Utils::getVector(char mode,
                                             Teuchos::RCP<Epetra_Vector> vec)
{
    if (mode == 'C') // copy
    {
        Teuchos::RCP<Epetra_Vector> copy =
            Teuchos::rcp(new Epetra_Vector(*vec));
        return copy;
    }
    else if (mode == 'V') // view
    {
        return vec;
    }
    else
    {
        WARNING("Invalid mode", __FILE__, __LINE__);
        return Teuchos::null;
    }
}

//=============================================================================
std::shared_ptr<std::vector<double> > Utils::getVector
(char mode, std::shared_ptr<std::vector<double> > vec)
{
    if (mode == 'C')      // copy
    {
        std::shared_ptr<std::vector<double> > copy =
            std::make_shared<std::vector<double> >(*vec);
        return copy;
    }
    else if (mode == 'V') // view
    {
        return vec;
    }
    else
    {
        WARNING("invalid mode", __FILE__, __LINE__);
        return std::shared_ptr<std::vector<double> >();
    }
}

//=============================================================================
void Utils::assembleCRS(Teuchos::RCP<Epetra_CrsMatrix> mat,
                        CRSMat const &crs, int const maxnnz,
                        Teuchos::RCP<TRIOS::Domain> domain)
{
    // we need domain information when the source crs is local
    bool global = (domain == Teuchos::null) ? true : false;

    // the beg array indicates whether the crs is 0 or 1-based
    bool index0;
    if (crs.beg[0] == 0)
        index0 = true;
    else if (crs.beg[0] == 1)
        index0 = false;
    else
    {
        WARNING("What CRS format is this? Continue with empty matrix.", __FILE__, __LINE__);
        return;
    }

    // indices array
    std::vector<int> indices(maxnnz, 0);

    // values array
    std::vector<double> values(maxnnz, 0.0);

    // define the rowmap we use
    Teuchos::RCP<Epetra_Map> rowMap;

    if (global)
    {
        rowMap = Teuchos::rcp(new Epetra_Map(mat->RowMap()));
        assert(rowMap->NumGlobalElements() == (int) crs.beg.size() - 1);
    }
    else
    {
        // in the case of a local crs assembly gives the GID's
        rowMap = domain->GetAssemblyMap();
        assert(rowMap->NumMyElements() == (int) crs.beg.size() - 1);
    }

#ifdef DEBUGGING_NEW
    // std::ofstream file;
    // file.open("rowmap" + std::to_string(rowMap->Comm().MyPID()));
    // rowMap->Print(file);
    // file.close();
#endif

    int numMyElements = rowMap->NumMyElements();

    int tRow, bRow, gRow, index, numEntries, col;
    int offset = (index0) ? 0 : 1;
    for (int lRow = 0; lRow < numMyElements; ++lRow)
    {
        // map using current map to global ID
        gRow = rowMap->GID(lRow);

        // map GID back through standardmap: if this gives -1 we have a
        // ghost point
        tRow = mat->RowMap().LID(gRow);

        if ( !global && (tRow == -1) )
        {
            continue;
        }

        bRow = (global) ? gRow : lRow;

        index      = crs.beg[bRow];
        numEntries = crs.beg[bRow+1] - index;

        // if we encounter a dense row (probably an integral equation)
        if (numEntries > maxnnz)
        {
            // adjust arrays
            indices = std::vector<int>(numEntries, 0);
            values  = std::vector<double>(numEntries, 0);
        }

        for (int j = 0; j < numEntries; ++j)
        {
            // taking 0-1-basedness into account using offset
            col = crs.jco[index + j - offset] - offset;

            indices[j] = (global) ? col : rowMap->GID(col);
            values[j]  = crs.co[index + j - offset];
        }

        int ierr;
        if (mat->Filled())
        {
            ierr =
                mat->ReplaceGlobalValues(gRow, numEntries,
                                         &values[0], &indices[0]);
        }
        else
        {
            ierr =
                mat->InsertGlobalValues(gRow, numEntries,
                                        &values[0], &indices[0]);
        }

        if (ierr != 0)
        {
            std::cout << "indices : ";
            for (int ii = 0; ii != numEntries; ++ii)
            {
                std::cout << indices[ii] << " ";
            }
            std::cout << std::endl;
            std::cout << "values : ";
            for (int ii = 0; ii != numEntries; ++ii)
            {
                std::cout << values[ii] << " ";
            }
            std::cout << std::endl;

            std::cout << "Error in Insert/ReplaceGlobalValues: "
                      << ierr << std::endl;

            std::cout << "Filled = " << mat->Filled()   << std::endl;
            std::cout << "Global = " << global          << std::endl;
            std::cout << "  GRID = " << gRow            << std::endl;
            std::cout << "  LRID = " << mat->LRID(gRow) << std::endl;

            std::cout << "jco : ";
            for (int jj = 0; jj != numEntries; ++jj)
            {
                col = crs.jco[index + jj - offset] - offset;
                std::cout << col << " ";
            }
            std::cout << std::endl;

            ERROR("Error in Insert/ReplaceGlobalValues",
                  __FILE__, __LINE__);
        }
    }
}

//=============================================================================
int Utils::SplitBox(int nx, int ny, int nz,
                    int nparts, int& ndx, int& ndy, int& ndz,
                    int sx, int sy, int sz)
{
    // Factor the number of processors into two dimensions. (nprocs = npN*npM)
    double rmin = 1e100;
    int ret = 1;

    int npx = nx / sx;
    int npy = ny / sy;
    int npz = nz / sz;

    std::string s1 = Teuchos::toString(nx) + "x" +
        Teuchos::toString(ny) + "x" + Teuchos::toString(nz);
    std::string s2 = Teuchos::toString(npx) + "x" +
        Teuchos::toString(npy) + "x" + Teuchos::toString(npz);

    // check all possibilities:
    for (int t1 = 1; t1 <= nparts; t1++)
        for (int t2 = 1; t2 <= (int)(nparts / t1); t2++)
        {
            int t3 = (int)(nparts/(t1*t2));

            if (t1 * t2 * t3 == nparts)
            {
                std::string s3 = Teuchos::toString(t1) + "x" +
                    Teuchos::toString(t2) + "x" + Teuchos::toString(t3);
                int my_nx = nx / t1;
                int my_ny = ny / t2;
                int my_nz = nz / t3;
                if ((my_nx * t1 != nx) || (my_ny * t2 != ny) || (my_nz * t3 != nz))
                {
                    DEBUG("Can't partition a "+s1+" domain into "+s3+" parts.");
                    continue;
                }
                int my_npx = npx / t1;
                int my_npy = npy / t2;
                int my_npz = npz / t3;
                if ((my_npx * sx != my_nx) ||
                    (my_npy * sy != my_ny) || (my_npz * sz != my_nz))
                {
                    DEBUG("Can't partition "+s2+" domains onto "+s3+" processors.");
                    continue;
                }

                double r1 = std::abs((double)nx / (double)t1 - (double)ny / (double)t2);
                double r2 = std::abs((double)nx / (double)t1 - (double)nz / (double)t3);
                double r3 = std::abs((double)ny / (double)t2 - (double)nz / (double)t3);
                double r = r1 + r2 + r3;

                if (r < rmin)
                {
                    rmin = r;
                    ndx  = t1;
                    ndy  = t2;
                    ndz  = t3;
                    ret  = 0;
                }
            }
        }
    return ret;
}
//========================================================================================
void Utils::ind2sub(int nx, int ny, int nz, int dof,
                    int idx, int& i, int& j, int& k, int& var)
{
#ifdef TESTING
    if ((idx < 0) || (idx >= (nx * ny * nz * dof)))
    {
        std::cerr << "dim=[" << nx << "," << ny << ","
                  << nz << "], dof=" << dof << ": ind=" <<idx << std::endl;
        ERROR("ind2sub: Index out of range!",__FILE__,__LINE__);
    }
#endif
    int rem = idx;
    var = MOD(rem, dof);
    rem = (rem - var) / dof;
    i   = MOD(rem , nx);
    rem = (rem - i) / nx;
    j   = MOD(rem, ny);
    rem = (rem - j) / ny;
    k   = MOD(rem, nz);
}
//========================================================================================
void Utils::ind2sub(int nx, int ny, int nz, int idx, int& i, int& j, int& k)
{
    int dummy;
    ind2sub(nx,ny,nz,1,idx,i,j,k,dummy);
    return;
}
//========================================================================================
//! converts cartesian subscripts to linear index
int Utils::sub2ind(int nx, int ny, int nz, int dof, int i, int j, int k, int var)
{
#ifdef TESTING
    std::string msg1 = "sub2ind: ";
    std::string msg3 = " out of range ";
    if ((i < 0) || (i >= nx))
    {
        std::string msg2 = "i-Index "+Teuchos::toString(i);
        std::string msg4 = "[0,"+Teuchos::toString(nx)+"]";
        ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
    }
    if ((j < 0) || ( j >= ny))
    {
        std::string msg2 = "j-Index "+Teuchos::toString(j);
        std::string msg4 = "[0,"+Teuchos::toString(ny)+"]";
        ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
    }
    if ((k < 0) || (k >= nz))
    {
        std::string msg2 = "k-Index "+Teuchos::toString(j);
        std::string msg4 = "[0,"+Teuchos::toString(nz)+"]";
        ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
    }
    if ((var<0)||(var>=dof))
    {
        std::string msg2 = "var-Index "+Teuchos::toString(j);
        std::string msg4 = "[0,"+Teuchos::toString(dof)+"]";
        ERROR(msg1+msg2+msg3+msg4,__FILE__,__LINE__);
    }
#endif
    return ((k*ny+j)*nx+i)*dof+var;
}

//========================================================================================
Teuchos::RCP<Epetra_Map> Utils::CreateMap(int nx, int ny, int nz,
                                          int dof, int indexbase,
                                          const Epetra_Comm& comm,
                                          int numActiveProcs)
{
    if (numActiveProcs == -1) numActiveProcs = comm.NumProc();
    // create a parallel map. We first figure out where in the domain we are
    int np  = std::min(comm.NumProc(),nx*ny*nz);
    np      = std::min(np,numActiveProcs);
    int pid = comm.MyPID();

    int npX, npY, npZ;
    int pidX  = -1, pidY = -1, pidZ = -1;
    int offX  = -1, offY = -1, offZ = -1;
    int nXloc = 0, nYloc = 0, nZloc = 0;

    bool active = (pid < np);

    Utils::SplitBox(nx, ny, nz, np, npX, npY, npZ);

    if (active)
    {
        Utils::ind2sub(npX, npY, npZ, pid, pidX, pidY, pidZ);
        // dimension of subdomain
        nXloc = (int) (nx / npX);
        nYloc = (int) (ny / npY);
        nZloc = (int) (nz / npZ);

        // offsets for local->global index conversion
        offX = pidX * (int) (nx / npX);
        offY = pidY * (int) (ny / npY);
        offZ = pidZ * (int) (nz / npZ);

        // distribute remaining points among first few cpu's
        int remX = nx % npX;
        int remY = ny % npY;
        int remZ = nz % npZ;

        if (pidX < remX) nXloc++;
        if (pidY < remY) nYloc++;
        if (pidZ < remZ) nZloc++;

        for (int i=0;i<std::min(remX,pidX);i++) offX++;
        for (int i=0;i<std::min(remY,pidY);i++) offY++;
        for (int i=0;i<std::min(remZ,pidZ);i++) offZ++;
    } //active?
    if (indexbase != 0)
    {
        // this could easily be implemented but I didn't need it up to now.
        ERROR("only index base 0 is implemented right now",
              __FILE__, __LINE__);
    }

    return CreateMap(offX, offX + nXloc - 1,
                     offY, offY + nYloc - 1,
                     offZ, offZ + nZloc - 1,
                     0, nx-1, 0, ny-1, 0, nz-1,
                     dof, comm);

}
//========================================================================================
Teuchos::RCP<Epetra_Map> Utils::CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
                                          int I0, int I1, int J0, int J1, int K0, int K1,
                                          int dof,
                                          const Epetra_Comm& comm)
{
    Teuchos::RCP<Epetra_Map> result = Teuchos::null;
    DEBUG("Utils::CreateMap ");
    DEBUG("["<<i0<<".."<<i1<<"]");
    DEBUG("["<<j0<<".."<<j1<<"]");
    DEBUG("["<<k0<<".."<<k1<<"]");

    int n = std::max(i1-i0+1,0); int N=I1-I0+1;
    int m = std::max(j1-j0+1,0); int M=J1-J0+1;
    int l = std::max(k1-k0+1,0); int L=K1-K0+1;

    DEBVAR(N);
    DEBVAR(M);
    DEBVAR(L);

    int NumMyElements = n*m*l*dof;
    int NumGlobalElements = -1; // note that there may be overlap
    int *MyGlobalElements = new int[NumMyElements];

    int pos = 0;
    for (int k=k0; k<=k1; k++)
        for (int j=j0; j<=j1; j++)
            for (int i=i0; i<=i1; i++)
                for (int var=0;var<dof;var++)
                {
                    MyGlobalElements[pos++] =
                        Utils::sub2ind(N,M,L,dof,i,j,k,var);
                }
    result = Teuchos::rcp(new Epetra_Map(NumGlobalElements,
                                         NumMyElements,MyGlobalElements,0,comm));
    delete [] MyGlobalElements;
    return result;
}
//========================================================================================
Teuchos::RCP<Epetra_Map> Utils::CreateMap(int i0, int i1, int j0, int j1, int k0, int k1,
                                          int I0, int I1, int J0, int J1, int K0, int K1,
                                          const Epetra_Comm& comm)
{
    Teuchos::RCP<Epetra_Map> result = Teuchos::null;

    DEBUG("Utils::CreateMap ");
    DEBUG("["<<i0<<".."<<i1<<"]");
    DEBUG("["<<j0<<".."<<j1<<"]");
    DEBUG("["<<k0<<".."<<k1<<"]");

    int n = i1-i0+1; int N=I1-I0+1;
    int m = j1-j0+1; int M=J1-J0+1;
    int l = k1-k0+1; //int L=K1-K0+1;

    DEBVAR(M);
    DEBVAR(N);
    DEBVAR(L);

    int NumMyElements = n*m*l;
    int NumGlobalElements = -1; // note that there may be overlap
    int *MyGlobalElements = new int[NumMyElements];

    int pos = 0;
    for (int k=k0; k<=k1; k++)
        for (int j=j0; j<=j1; j++)
            for (int i=i0; i<=i1; i++)
                MyGlobalElements[pos++] = k*N*M + j*N + MOD((double)i,(double)N);

    result = Teuchos::rcp(new Epetra_Map(NumGlobalElements,
                                         NumMyElements,MyGlobalElements,0,comm));
    delete [] MyGlobalElements;
    return result;
}
//========================================================================================
//! extract a map with nun = nvars from a map with nun=6. 'var'
//! is the array of variables to be extracted.
Teuchos::RCP<Epetra_Map> Utils::CreateSubMap(const Epetra_Map& map,
                                             int dof, int var)
{
    return CreateSubMap(map,dof,&var,1);
}
//========================================================================================
//! extract a map with nun=2 from a map with nun=6. 'var'
//! are the variables to be extracted, i.e. {UU,VV}, {TT,SS} etc.
Teuchos::RCP<Epetra_Map> Utils::CreateSubMap(const Epetra_Map& map,
                                             int dof, const int var[2])
{
    return CreateSubMap(map,dof,var,2);
}

//========================================================================================
//! extract a map with nun = nvars from a map with nun = 6. 'var'
//! is the array of variables to be extracted.
Teuchos::RCP<Epetra_Map> Utils::CreateSubMap(const Epetra_Map& map,
                                             int dof, const int *var, int nvars)
{
    int dim    = map.NumMyElements();     // number of entries in original map
    int dimGlb = map.NumGlobalElements(); // number of global entries in original map

    int numBlocks     = dim/dof;          // number of blocks
    int numBlocksGlb  = dimGlb/dof;       // global number of blocks


    if (numBlocks * dof < dim-1)
    {
        ERROR("\nInvalid dimension detected, possibly more \n" <<
              " than one auxiliary equation in map... ",
              __FILE__, __LINE__);
    }


    // Handle possible auxiliary variables
    bool auxVar[nvars];
    int  auxRws[nvars];
    for (int j = 0; j < nvars; j++)
    {
        if (var[j] > dof)
        {
            auxVar[j] = true;
            auxRws[j] = numBlocksGlb * dof + (var[j] - dof - 1);

            std::cout << " pid " << map.Comm().MyPID()
                      << " auxiliary variable detected: " << var[j]
                      << ", row " << auxRws[j] << std::endl;
        }
        else
        {
            auxVar[j] = false;
            auxRws[j] = -1;
        }
    }

    int subdim = numBlocks * nvars;   // number of local entries in new map (<=dim)
    int *MyGlobalElements = new int[subdim];

    // take the entries from the old map that correspond
    // to those in 'vars' and put them in the input array
    // for the new map.

    int k = 0;
    for (int i  = 0; i < numBlocks; i++)
        for (int j = 0; j < nvars; j++)
        {
            if (!auxVar[j])
                MyGlobalElements[k++] = map.GID(i*dof+(var[j]-1));
        }

    for (int j = 0; j < nvars; j++)
    {
        if (auxVar[j])
            MyGlobalElements[k++] = auxRws[j];
    }

    Teuchos::RCP<Epetra_Map> submap =
        Teuchos::rcp(new Epetra_Map(-1, k, MyGlobalElements, 0, map.Comm()));

    delete [] MyGlobalElements;
    return submap;
}

//========================================================================================
//! given a map and an array indicating wether each node of the map is to be
//! discarded (true) or not (false), this function creates a new map with the
//! discarded entries removed.
Teuchos::RCP<Epetra_Map> Utils::CreateSubMap(const Epetra_Map& map,
                                             const bool* discard)
{
    int numel = map.NumMyElements();
    int *MyGlobalElements = new int[numel]; // 'worst' case: no discarded nodes
    int numel_new = 0;
    for (int k=0;k<numel;k++)
    {
        if (!discard[k])
        {
            MyGlobalElements[numel_new] = map.GID(k);
            numel_new++;
        }
    }
    Teuchos::RCP<Epetra_Map> submap =
        Teuchos::rcp(new Epetra_Map(-1, numel_new, MyGlobalElements,
                                    map.IndexBase(), map.Comm()));
    delete [] MyGlobalElements;
    return submap;
}


//========================================================================================
//! Given a map and a list of global indices we create a submap with the same parallel
//! distribution but restricted to the given indices.
//! --> this has not been properly tested yet!!
Teuchos::RCP<Epetra_BlockMap> Utils::CreateSubMap(const Epetra_BlockMap& map,
                                                  std::vector<int> const &list)
{
    int listdim = list.size();             // number of entries in new map

    if (listdim > map.NumGlobalElements()) // weak check of correct list
        ERROR("unexpected number of elements in index list!",__FILE__,__LINE__);

    std::vector<int> MyGlobalElements;  // Container for elements owned

    int *PIDList = new int[listdim];    // List of PID's
    int *LIDList = new int[listdim];    // List of LID's

    // Get a list of PID's and LID's in the map based on the list of global
    // indices from the given arguments.
    map.RemoteIDList(listdim, &list[0], PIDList, LIDList);

    // Create MyGlobalElements array based on the information in PIDList
    int myPID = map.Comm().MyPID();
    for (int i = 0; i < listdim; ++i)
        if (PIDList[i] == myPID)
            MyGlobalElements.push_back(list[i]);

    // Create submap
    Teuchos::RCP<Epetra_Map> submap =
        Teuchos::rcp(new Epetra_Map(listdim, MyGlobalElements.size(),
                                    &MyGlobalElements[0], 0, map.Comm()));

    delete [] PIDList;
    delete [] LIDList;
    return submap;
}

//=======================================================================================
//! Given an Epetra_Vector and a list of global indices we return an Epetra_Vector with
//! the same distribution but restricted to the values at the supplied indices.
//! ITS SLOW!
Teuchos::RCP<Epetra_Vector> Utils::RestrictVector(Epetra_Vector const &vector,
                                                  std::vector<int> const &indices)
{
    TIMER_START("Utils: restrict vector");
    // Create restricted map
    Teuchos::RCP<Epetra_BlockMap> indexMap =
        Utils::CreateSubMap(vector.Map(), indices);

    // Create the output vector
    Teuchos::RCP<Epetra_Vector> restrictedVector =
        Teuchos::rcp(new Epetra_Vector(*indexMap) );

    TIMER_STOP("Utils: restrict vector");

    TIMER_START("Utils: restrict vector import");

    // Create importer
    // target map: indexMap
    // source map: vector.Map()
    Epetra_Import full2restricted(*indexMap, vector.Map());

    // Import vector values into restricted vector
    restrictedVector->Import(vector, full2restricted, Insert);
    TIMER_STOP("Utils: restrict vector import");

    return restrictedVector;
}

//========================================================================================
Teuchos::RCP<Epetra_CrsMatrix> Utils::Gather(const Epetra_CrsMatrix& mat, int root)
{
    const Epetra_Map& rowmap_dist = mat.RowMap();
    // we take the domain map as the colmap is potentially overlapping
    const Epetra_Map& colmap_dist = mat.DomainMap();
    // gather the row map
    Teuchos::RCP<Epetra_Map> rowmap = Gather(rowmap_dist,root);
    // gather the col map
    Teuchos::RCP<Epetra_Map> colmap = Gather(colmap_dist,root);

    // we only guess the number of row entries, this routine is not performance critical
    // as it should only be used for debugging anyway
    int num_entries = mat.NumGlobalNonzeros() / mat.NumGlobalRows();
    Teuchos::RCP<Epetra_CrsMatrix> gmat =
        Teuchos::rcp(new Epetra_CrsMatrix(Copy,*rowmap, *colmap, num_entries) );

    Teuchos::RCP<Epetra_Import> import =
        Teuchos::rcp(new Epetra_Import(*rowmap,rowmap_dist) );
    CHECK_ZERO(gmat->Import(mat,*import,Insert));
    CHECK_ZERO(gmat->FillComplete());
    gmat->SetLabel(mat.Label());
    return gmat;
}

//========================================================================================
Teuchos::RCP<Epetra_MultiVector> Utils::Gather(const Epetra_MultiVector& vec, int root)
{

    const Epetra_BlockMap& map_dist = vec.Map();
    Teuchos::RCP<Epetra_BlockMap> map = Gather(map_dist,root);
    Teuchos::RCP<Epetra_MultiVector> gvec =
        Teuchos::rcp(new Epetra_MultiVector(*map,vec.NumVectors()));
    Teuchos::RCP<Epetra_Import> import =
        Teuchos::rcp(new Epetra_Import(*map,map_dist) );

    CHECK_ZERO(gvec->Import(vec,*import,Insert));
    gvec->SetLabel(vec.Label());
    return gvec;
}

//========================================================================================
// create "Gather" map from "Solve" map
Teuchos::RCP<Epetra_BlockMap> Utils::Gather(const Epetra_BlockMap& map, int root)
{

    int ElementSize = map.ElementSize();
#ifdef TESTING
    if (ElementSize != 1)
    {
        DEBVAR(ElementSize);
        ERROR("this is possibly not implemented correctly!",
              __FILE__, __LINE__);
        ElementSize=1;
    }
#endif
    int NumMyElements       = map.NumMyElements();
    int NumGlobalElements   = map.NumGlobalElements();
    const Epetra_Comm& Comm = map.Comm();
    int *MyGlobalElements   = new int[NumMyElements];
    int *AllGlobalElements  = NULL;


    for (int i = 0; i < NumMyElements; i++)
    {
        MyGlobalElements[i] = map.GID(i);
    }
    if (Comm.MyPID() == root)
    {
        AllGlobalElements = new int[NumGlobalElements];
    }
    if (Comm.NumProc()>1)
    {
#ifdef HAVE_MPI
        const Epetra_MpiComm MpiComm = dynamic_cast<const Epetra_MpiComm&>(Comm);
        int *counts, *disps;
        counts = new int[Comm.NumProc()];
        disps  = new int[Comm.NumProc()+1];
        MPI_Gather(&NumMyElements,1,MPI_INTEGER,
                   counts,1,MPI_INTEGER,root,MpiComm.GetMpiComm());

        if (Comm.MyPID()==root)
        {
            disps[0]=0;
            for (int p=0;p<Comm.NumProc();p++)
            {
                disps[p+1] = disps[p]+counts[p];
            }
        }

        MPI_Gatherv(MyGlobalElements, NumMyElements,MPI_INTEGER,
                    AllGlobalElements, counts,disps, MPI_INTEGER, root, MpiComm.GetMpiComm());
        delete [] counts;
        delete [] disps;
#else
        ERROR("No MPI but still parallel??? We don't do that.", __FILE__, __LINE__);
#endif
    }
    else
    {
        for (int i = 0; i < NumMyElements; i++)
            AllGlobalElements[i] = MyGlobalElements[i];
    }

    if (Comm.MyPID()!=root)
    {
        NumMyElements=0;
    }
    else
    {
        NumMyElements = NumGlobalElements;
        Teuchos::ArrayView<int> view(AllGlobalElements, NumGlobalElements);
        std::sort(view.begin(), view.end());
        Teuchos::ArrayView<int>::iterator new_end = std::unique(view.begin(), view.end());
        NumMyElements = std::distance(view.begin(), new_end);
        NumGlobalElements = NumMyElements;
    }
    CHECK_ZERO(Comm.Broadcast(&NumGlobalElements, 1, root));
    // build the new (gathered) map
    Teuchos::RCP<Epetra_BlockMap> gmap =
        Teuchos::rcp(new Epetra_BlockMap(NumGlobalElements, NumMyElements, AllGlobalElements,
                                         ElementSize, map.IndexBase(), Comm) );


    if (Comm.MyPID() == root)
    {
        delete [] AllGlobalElements;
    }
    delete [] MyGlobalElements;

    return gmap;
}
//========================================================================================
// create "Gather" map from "Solve" map
Teuchos::RCP<Epetra_Map> Utils::Gather(const Epetra_Map& map, int root)
{
    int NumMyElements = map.NumMyElements();
    int NumGlobalElements = map.NumGlobalElements();
    const Epetra_Comm& Comm = map.Comm();

    int *MyGlobalElements = new int[NumMyElements];
    int *AllGlobalElements = NULL;

    for (int i = 0; i < NumMyElements; i++)
        MyGlobalElements[i] = map.GID(i);

    if (Comm.MyPID() == root)
        AllGlobalElements = new int[NumGlobalElements];

    if (Comm.NumProc()>1)
    {
#ifdef HAVE_MPI
        const Epetra_MpiComm MpiComm = dynamic_cast<const Epetra_MpiComm&>(Comm);
        int *counts, *disps;
        counts = new int[Comm.NumProc()];
        disps  = new int[Comm.NumProc()+1];
        MPI_Gather(&NumMyElements, 1, MPI_INTEGER,
                   counts, 1, MPI_INTEGER, root, MpiComm.GetMpiComm());

        if (Comm.MyPID() == root)
        {
            disps[0] = 0;
            for (int p = 0; p < Comm.NumProc(); p++)
                disps[p+1] = disps[p] + counts[p];
        }

        MPI_Gatherv(MyGlobalElements,  NumMyElements,MPI_INTEGER,
                    AllGlobalElements, counts,disps, MPI_INTEGER, root, MpiComm.GetMpiComm());
        delete [] counts;
        delete [] disps;
#else
        ERROR("No MPI but still parallel??? We don't do that.",__FILE__,__LINE__);
#endif
    }
    else
        for (int i=0; i < NumMyElements; i++)
            AllGlobalElements[i] = MyGlobalElements[i];

    if (Comm.MyPID() != root)
        NumMyElements = 0;
    else
    {
        NumMyElements=NumGlobalElements;
        std::sort(AllGlobalElements,AllGlobalElements+NumGlobalElements);
    }
    // build the new (gathered) map
    Teuchos::RCP<Epetra_Map> gmap =
        Teuchos::rcp(new Epetra_Map(NumGlobalElements, NumMyElements, AllGlobalElements,
                                    map.IndexBase(), Comm));

    if (Comm.MyPID() == root)
        delete [] AllGlobalElements;

    delete [] MyGlobalElements;

    return gmap;
}
//========================================================================================
// distribute a gathered vector among processors
Teuchos::RCP<Epetra_MultiVector> Utils::Scatter
(const Epetra_MultiVector& vec, const Epetra_BlockMap& distmap)
{
    Teuchos::RCP<Epetra_MultiVector> dist_vec =
        Teuchos::rcp(new Epetra_MultiVector(distmap,vec.NumVectors()));
    Teuchos::RCP<Epetra_Import> import =
        Teuchos::rcp(new Epetra_Import(vec.Map(),distmap));
    CHECK_ZERO(dist_vec->Export(vec,*import,Insert));
    return dist_vec;
}
//========================================================================================
// create "col" map from "Solve" map
Teuchos::RCP<Epetra_BlockMap> Utils::AllGather(const Epetra_BlockMap& map, bool reorder)
{
    int ElementSize         = map.ElementSize();
    int NumMyElements       = map.NumMyElements();
    int NumGlobalElements   = map.NumGlobalElements();
    const Epetra_Comm& Comm = map.Comm();
#ifdef TESTING
    if (ElementSize > 1)
    {
        ERROR("this is possibly not implemented correctly!",
              __FILE__, __LINE__);
    }
#endif
    int *MyGlobalElements  = new int[NumMyElements];
    int *AllGlobalElements = new int[NumGlobalElements];
    for (int i = 0; i < NumMyElements; ++i)
    {
        MyGlobalElements[i] = map.GID(i);
    }
    if (Comm.NumProc() > 1)
    {
#ifdef HAVE_MPI
        const Epetra_MpiComm MpiComm = dynamic_cast<const Epetra_MpiComm&>(Comm);
        int *counts, *disps;
        counts = new int[Comm.NumProc()];
        disps = new int[Comm.NumProc()+1];
        MPI_Allgather(&NumMyElements,1,MPI_INTEGER,
                      counts,1,MPI_INTEGER,MpiComm.GetMpiComm());

        disps[0]=0;
        for (int p=0;p<Comm.NumProc();p++)
        {
            disps[p+1] = disps[p]+counts[p];
        }

        MPI_Allgatherv(MyGlobalElements, NumMyElements,MPI_INTEGER,
                       AllGlobalElements, counts,disps, MPI_INTEGER, MpiComm.GetMpiComm());
        delete [] counts;
        delete [] disps;
#else
        ERROR("No MPI but still parallel? We don't do tthat.",__FILE__,__LINE__);
#endif
    }
    else
    {
        for (int i = 0; i < NumMyElements; ++i)
            AllGlobalElements[i] = MyGlobalElements[i];
    }

    NumMyElements = NumGlobalElements;
    NumGlobalElements = -1;

    if (reorder)
    {
        std::sort(AllGlobalElements,AllGlobalElements+NumMyElements);
    }
    // build the new (gathered) map
    Teuchos::RCP<Epetra_BlockMap> gmap =
        Teuchos::rcp(new Epetra_BlockMap (NumGlobalElements, NumMyElements,
                                          AllGlobalElements, ElementSize, map.IndexBase(), Comm) );

    delete [] MyGlobalElements;
    delete [] AllGlobalElements;

    return gmap;
} //AllGather 1
//========================================================================================
// create "col" map from "Solve" map
Teuchos::RCP<Epetra_Map> Utils::AllGather(const Epetra_Map& map, bool reorder)
{
    int NumMyElements       = map.NumMyElements();
    int NumGlobalElements   = map.NumGlobalElements();
    const Epetra_Comm& Comm = map.Comm();

    int *MyGlobalElements  = new int[NumMyElements];
    int *AllGlobalElements = new int[NumGlobalElements];

    for (int i = 0; i < NumMyElements; ++i)
    {
        MyGlobalElements[i] = map.GID(i);
    }

    if (Comm.NumProc() > 1)
    {
#ifdef HAVE_MPI
        const Epetra_MpiComm MpiComm = dynamic_cast<const Epetra_MpiComm&>(Comm);
        int *counts, *disps;
        counts = new int[Comm.NumProc()];
        disps  = new int[Comm.NumProc()+1];
        MPI_Allgather(&NumMyElements,1,MPI_INTEGER,
                      counts,1,MPI_INTEGER,MpiComm.GetMpiComm());

        disps[0] = 0;
        for (int p = 0; p < Comm.NumProc(); ++p)
        {
            disps[p+1] = disps[p] + counts[p];
        }

        MPI_Allgatherv(MyGlobalElements, NumMyElements,MPI_INTEGER,
                       AllGlobalElements, counts,disps, MPI_INTEGER, MpiComm.GetMpiComm());
        delete [] counts;
        delete [] disps;
#else
        ERROR("No MPI but still parallel? We don't do tthat.",__FILE__,__LINE__);
#endif
    }
    else
    {
        for (int i = 0; i < NumMyElements; i++)
            AllGlobalElements[i] = MyGlobalElements[i];
    }

    NumMyElements = NumGlobalElements;
    NumGlobalElements = -1;

    if (reorder)
    {
        std::sort(AllGlobalElements,AllGlobalElements+NumMyElements);
    }

    // build the new (gathered) map
    Teuchos::RCP<Epetra_Map> gmap =
        Teuchos::rcp(new Epetra_Map (NumGlobalElements, NumMyElements,
                                     AllGlobalElements, map.IndexBase(), Comm) );

    delete [] MyGlobalElements;
    delete [] AllGlobalElements;

    return gmap;
} //AllGather
//========================================================================================
Teuchos::RCP<Epetra_MultiVector> Utils::AllGather(const Epetra_MultiVector& vec)
{
    TIMER_START("Utils: all gather");
    const Epetra_BlockMap& map_dist = vec.Map();
    Teuchos::RCP<Epetra_BlockMap> map = AllGather(map_dist);
    Teuchos::RCP<Epetra_MultiVector> gvec = Teuchos::rcp(new Epetra_Vector(*map,vec.NumVectors()));
    TIMER_STOP("Utils: all gather");
    TIMER_START("Utils: all gather create import");
    Teuchos::RCP<Epetra_Import> import = Teuchos::rcp(new Epetra_Import(*map,map_dist) );
    TIMER_STOP("Utils: all gather create import");
    TIMER_START("Utils: all gather import");
    gvec->Import(vec,*import,Insert);
    TIMER_STOP("Utils: all gather import");
    gvec->SetLabel(vec.Label());
    return gvec;
}
//========================================================================================
Teuchos::RCP<Epetra_IntVector> Utils::AllGather(const Epetra_IntVector& vec)
{
    const Epetra_BlockMap& map_dist = vec.Map();
    Teuchos::RCP<Epetra_BlockMap> map = AllGather(map_dist);
    Teuchos::RCP<Epetra_IntVector> gvec = Teuchos::rcp(new Epetra_IntVector(*map));
    Teuchos::RCP<Epetra_Import> import = Teuchos::rcp(new Epetra_Import(*map,map_dist) );
    gvec->Import(vec,*import,Insert);
    gvec->SetLabel(vec.Label());
    return gvec;
}
//========================================================================================
Teuchos::RCP<Epetra_CrsMatrix> Utils::MatrixProduct(bool transA, const Epetra_CrsMatrix& A,
                                                    bool transB, const Epetra_CrsMatrix& B,
                                                    bool useColMap)
{
    DEBUG("Entered MatrixProduct()");

    DEBUG("  construct AB");
    Teuchos::RCP<Epetra_CrsMatrix> AB;
    DEBVAR(useColMap);
    if (useColMap)
    {
        AB = Teuchos::rcp(new Epetra_CrsMatrix(Copy, A.RowMap(),
                                               B.ColMap(), A.MaxNumEntries()));
        DEBVAR(B.Importer());
        DEBVAR(B.ColMap().SameAs(AB->ColMap()));
    }
    else
    {
        AB = Teuchos::rcp(new Epetra_CrsMatrix(Copy, A.RowMap(),
                                               A.MaxNumEntries()));
    }

    DEBUG("  compute A*B...");
    DEBVAR(transA);
    DEBVAR(A.NumGlobalRows());
    DEBVAR(A.NumGlobalCols());
    DEBVAR(transB);
    DEBVAR(B.NumGlobalRows());
    DEBVAR(B.NumGlobalCols());

#ifdef TESTING
    if (!A.Filled())
    {
        ERROR("Matrix A not filled!",__FILE__,__LINE__);
    }
    else
    {
        DEBUG("Matrix A filled!");
    }
    if (!B.Filled())
    {
        ERROR("Matrix B not filled!",__FILE__,__LINE__);
    }
    else
    {
        DEBUG("Matrix B filled!");
    }
#endif

    DEBUG(" perform A*B");
    DEBUG("--barrier sync");
    A.Comm().Barrier();
    CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(A,transA,B,transB,*AB));
    DEBUG("done!");
    DEBUG("Leaving MatrixProduct()");
    return AB;
}
//========================================================================================
// compress a matrix' column map so that the resulting map contains
// only points actually appearing as column indices of the matrix
Teuchos::RCP<Epetra_Map> Utils::CompressColMap(const Epetra_CrsMatrix& A)
{
    if (!A.HaveColMap()) ERROR("Matrix has no column map!",__FILE__,__LINE__);
    const Epetra_Map& old_map = A.ColMap();
    int n_old = old_map.NumMyElements();
    bool *is_col_entry = new bool[n_old];

    for (int i = 0; i < n_old; i++)
        is_col_entry[i] = false;

    for (int i = 0; i < A.NumMyRows(); i++)
    {
        int *ind;
        int len;
        CHECK_ZERO(A.Graph().ExtractMyRowView(i,len,ind));
        for (int j=0;j<len;j++)
            is_col_entry[ind[j]] = true;
    }

    int n_new = 0;
    int *new_elements = new int[n_old];

    for (int i = 0; i < n_old; i++)
    {
        if (is_col_entry[i])
        {
            new_elements[n_new++] = old_map.GID(i);
        }
    }
    Teuchos::RCP<Epetra_Map> new_map =
        Teuchos::rcp(new Epetra_Map(-1,n_new,new_elements,old_map.IndexBase(),old_map.Comm()));

    delete [] new_elements;
    delete [] is_col_entry;

    return new_map;
}
//========================================================================================
// create an exact copy of a matrix replacing the column map.
// The column maps have to be 'compatible'
// in the sense that the new ColMap is a subset of the old one.
Teuchos::RCP<Epetra_CrsMatrix> Utils::ReplaceColMap(Teuchos::RCP<Epetra_CrsMatrix> A,
                                                    const Epetra_Map& newcolmap)
{
    int maxlen = A->MaxNumEntries();
    int len;
    int *ind = new int[maxlen];
    double *val = new double[maxlen];
    int nloc = A->NumMyRows();
    int *row_lengths = new int[nloc];
    for (int i=0;i<nloc;i++) row_lengths[i]=A->NumMyEntries(i);
    Teuchos::RCP<Epetra_CrsMatrix> tmpmat;
    tmpmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RowMap(),
                                               newcolmap, row_lengths) );
    int grid;
    for (int i=0;i<nloc;i++)
    {
        grid = A->GRID(i);
        CHECK_ZERO(A->ExtractGlobalRowCopy(grid,maxlen,len,val,ind));
#ifdef DEBUGGING
//      (*debug) << "row " << grid << ": ";
//      for (int j=0;j<len;j++) (*debug) << ind[j] << " ";
//      (*debug) << std::endl;
#endif
        CHECK_ZERO(tmpmat->InsertGlobalValues(grid, len, val, ind));
    }
    tmpmat->SetLabel(A->Label());
    delete [] ind;
    delete [] val;
    delete [] row_lengths;
    return tmpmat;
}
//========================================================================================
// workaround for the buggy Trilinos routine with the same name
Teuchos::RCP<Epetra_CrsMatrix> Utils::ReplaceRowMap(Teuchos::RCP<Epetra_CrsMatrix> A,
                                                    const Epetra_Map& newmap)
{
    int maxlen       = A->MaxNumEntries();
    int    *ind      = new int[maxlen];
    double *val      = new double[maxlen];
    int  nloc        = A->NumMyRows();
    int *row_lengths = new int[nloc];

    for (int i=0; i < nloc; i++)
        row_lengths[i] = A->NumMyEntries(i);

    Teuchos::RCP<Epetra_CrsMatrix> tmpmat;
    if (A->HaveColMap())
    {
        tmpmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,newmap,A->ColMap(), row_lengths) );
    }
    else
    {
        tmpmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,newmap, row_lengths) );
    }

    int len;
    int rowA,rowNew;
    for (int i=0;i<A->NumMyRows();i++)
    {
        rowA = A->GRID(i);
        rowNew = newmap.GID(i);
        CHECK_ZERO(A->ExtractGlobalRowCopy(rowA,maxlen,len,val,ind));
        CHECK_ZERO(tmpmat->InsertGlobalValues(rowNew, len, val, ind));
    }
    tmpmat->SetLabel(A->Label());
    delete [] ind;
    delete [] val;
    delete [] row_lengths;
    return tmpmat;
}
//========================================================================================
// create an exact copy of a matrix removing the column map.
// This means that row- and column map have to be 'compatible'
// in the sense that the ColMap is a subset of the RowMap.
// It seems to be required in order to use Ifpack in some cases.
Teuchos::RCP<Epetra_CrsMatrix> Utils::RemoveColMap(Teuchos::RCP<Epetra_CrsMatrix> A)
{
    int maxlen       = A->MaxNumEntries();
    int *ind         = new int[maxlen];
    double *val      = new double[maxlen];
    int nloc         = A->NumMyRows();
    int *row_lengths = new int[nloc];

    for (int i = 0; i < nloc; i++)
        row_lengths[i] = A->NumMyEntries(i);
    Teuchos::RCP<Epetra_CrsMatrix> tmpmat;
    tmpmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RowMap(), row_lengths));

    int len;
    int grid;
    for (int i = 0; i < A->NumMyRows(); i++)
    {
        grid = A->GRID(i);
        CHECK_ZERO(A->ExtractGlobalRowCopy(grid,maxlen,len,val,ind));
        CHECK_ZERO(tmpmat->InsertGlobalValues(grid, len, val, ind));
    }
    tmpmat->SetLabel(A->Label());
    delete [] ind;
    delete [] val;
    delete [] row_lengths;
    return tmpmat;
}
//========================================================================================
Teuchos::RCP<Epetra_CrsMatrix> Utils::RebuildMatrix(Teuchos::RCP<Epetra_CrsMatrix> A)
{
    int maxlen       = A->MaxNumEntries();
    int *ind         = new int[maxlen];
    double *val      = new double[maxlen];
    int nloc         = A->NumMyRows();
    int *row_lengths = new int[nloc];

    for (int i = 0; i < nloc; i++)
        row_lengths[i] = A->NumMyEntries(i);
    Teuchos::RCP<Epetra_CrsMatrix> tmpmat;
    tmpmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RowMap(), A->ColMap(), row_lengths));
    int len;
    int grid;
    for (int i = 0; i < A->NumMyRows(); i++)
    {
        grid = A->GRID(i);
        CHECK_ZERO(A->ExtractGlobalRowCopy(grid,maxlen,len,val,ind));
        CHECK_ZERO(tmpmat->InsertGlobalValues(grid, len, val, ind));
    }
    tmpmat->SetLabel(A->Label());
    delete [] ind;
    delete [] val;
    delete [] row_lengths;
    return tmpmat;
}
//========================================================================================
// simultaneously replace row and column map
Teuchos::RCP<Epetra_CrsMatrix> Utils::ReplaceBothMaps(Teuchos::RCP<Epetra_CrsMatrix> A,
                                                      const Epetra_Map& newmap,
                                                      const Epetra_Map& newcolmap)
{
    //DEBVAR(A->RowMap());
    //DEBVAR(newmap);
    //DEBVAR(A->ColMap());
    //DEBVAR(newcolmap);
    int maxlen = A->MaxNumEntries();
    int len;
    int    *ind = new int[maxlen];
    double *val = new double[maxlen];
    int nloc = A->NumMyRows();
    int *row_lengths = new int[nloc];
    for (int i=0;i<nloc;i++)
        row_lengths[i]=A->NumMyEntries(i);
    Teuchos::RCP<Epetra_CrsMatrix> tmpmat;
    tmpmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,newmap,newcolmap,row_lengths) );
    int rowA,rowNew;

    for (int i = 0; i < A->NumMyRows(); i++)
    {
        rowA = A->GRID(i);
        rowNew = newmap.GID(i);
        CHECK_ZERO(A->ExtractGlobalRowCopy(rowA,maxlen,len,val,ind));
        for (int j=0;j<len;j++)
        {
            int newind=newcolmap.GID(A->LCID(ind[j]));
//        DEBUG(i<<" ("<<rowA<<"->"<<rowNew<<"), "<<A->LCID(ind[j])<<"("<<ind[j]<<"->"<<newind<<")");
            ind[j] = newind;
        }
        CHECK_ZERO(tmpmat->InsertGlobalValues(rowNew, len, val, ind));
    }
    tmpmat->SetLabel(A->Label());
    delete [] ind;
    delete [] val;
    delete [] row_lengths;
    return tmpmat;
}
//========================================================================================
Teuchos::RCP<Epetra_CrsMatrix> Utils::TripleProduct(bool transA, const Epetra_CrsMatrix& A,
                                                    bool transB, const Epetra_CrsMatrix& B,
                                                    bool transC, const Epetra_CrsMatrix& C)
{
    // trans(A) is not available as we prescribe the row-map of A*B, but if it is needed
    // at some point it can be readily implemented
    if(transA) ERROR("This case is not implemented: trans(A)*op(B)*op(C)\n",__FILE__,__LINE__);

    // temp matrix
    Teuchos::RCP<Epetra_CrsMatrix> AB =
        Teuchos::rcp(new Epetra_CrsMatrix(Copy,A.RowMap(),A.MaxNumEntries()) );

    DEBUG("compute A*B...");
    CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(A,transA,B,transB,*AB));

    // result matrix
    Teuchos::RCP<Epetra_CrsMatrix> ABC =
        Teuchos::rcp(new Epetra_CrsMatrix(Copy,AB->RowMap(),AB->MaxNumEntries()) );

    DEBUG("compute ABC...");
    CHECK_ZERO(EpetraExt::MatrixMatrix::Multiply(*AB,false,C,transC,*ABC));

    DEBUG("done!");
    return ABC;
}
//========================================================================================
Teuchos::RCP<Epetra_Map> Utils::ExtractRange(const Epetra_Map& M, int i1, int i2)
{

#ifdef TESTING
    int n = M.MaxAllGID();
    if (i1<0||i1>n) ERROR("CreateSubMap: lower bound out of range!",__FILE__,__LINE__);
    if (i2<0||i2>n) ERROR("CreateSubMap: upper bound out of range!",__FILE__,__LINE__);
    if (i2<i1)      ERROR("CreateSubMap: invalid interval bounds!" ,__FILE__,__LINE__);
#endif

    int *MyGlobalElements = new int[M.NumMyElements()];
    int p=0;
    int gid;
    for (int i=0;i<M.NumMyElements();i++)
    {
        gid = M.GID(i);
        if (gid>=i1 && gid<=i2) MyGlobalElements[p++]=gid;
    }
    // build the two new maps. Set global num el. to -1 so Epetra recomputes it
    Teuchos::RCP<Epetra_Map> M1 =
        Teuchos::rcp(new Epetra_Map(-1,p,MyGlobalElements,M.IndexBase(),M.Comm()) );
    delete [] MyGlobalElements;
    return M1;
}
//========================================================================================
