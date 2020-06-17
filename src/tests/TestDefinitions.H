#ifndef TESTDEFINITIONS_H
#define TESTDEFINITIONS_H

//=============================================================================
// General definitions for tests
//=============================================================================

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "gtest/gtest.h"

#include <functional>
#include <vector>
#include <memory>

class Epetra_Vector;

// We need this to compute norms
extern "C" double ddot_(int *N, double *X, int *INCX, double *Y, int *INCY);

//------------------------------------------------------------------
double norm(Teuchos::RCP<Epetra_Vector> vec);

double norm(std::shared_ptr<std::vector<double> > vec);

std::shared_ptr<std::vector<double> >
getGatheredVector(Teuchos::RCP<Epetra_Vector> vec);

//------------------------------------------------------------------
// test the entries of numerical jacobian
template<typename ModelPtr, typename VectorPtr, typename CCS>
void testEntries(ModelPtr model, CCS &ccs, VectorPtr x)
{
    x->PutScalar(0.0);
    int myLength = x->MyLength();
    EXPECT_EQ(myLength, (int) (ccs.beg.size() - 1));
            
    typename VectorPtr::element_type e_i(*x);
    typename VectorPtr::element_type e_j(*x);
    typename VectorPtr::element_type tmp(*x);
            
    // assuming we run on a single core
    int    ico  = 0;
    double Jij  = 0.0;

    double error, diff;
    for (int jco = 0; jco != myLength; ++jco)
    {
        e_j[jco] = 1;
        tmp.PutScalar(0.0);

        // tmp = jth column 
        model->applyMatrix(e_j, tmp);
        for (int i = ccs.beg[jco]; i != ccs.beg[jco+1]; ++i)
        {
            ico = ccs.ico[i];
            e_i[ico] = 1;

            // get Jij component
            e_i.Dot(tmp, &Jij);
            e_i[ico] = 0;

            // difference between analytical and numerical value
            diff  = std::abs(Jij - ccs.co[i]);                    
                    
            error = std::max(1e-2 * std::abs(Jij), 1e-10);
            EXPECT_NEAR(diff, 0.0, error);
                                                  
            if (diff > error) // some additional outputting
            {
                std::cout << "(" << ico << "," << jco << ") ";
                std::cout << "nuJac: " << ccs.co[i] << " " << "anJac: " << Jij << std::endl;
                int mdl,i,j,k,xx;
                        
                model->gid2coord(ico,mdl,i,j,k,xx);
                std::cout << "\nrow: mdl=" << mdl << " i=" << i << " j=" << j
                          << " k=" << k << " xx=" << xx << std::endl;
                        
                model->gid2coord(jco,mdl,i,j,k,xx);
                std::cout << "col: mdl=" << mdl << " i=" << i << " j=" << j
                          << " k=" << k << " xx=" << xx << std::endl << std::endl;
            }
        }
        e_j[jco] = 0; // reset value
    }            
}

//------------------------------------------------------------------
// Functions for checking Teuchos::ParameterList's


// Combines two Google Test AssertionResults together by combining the messages
// of failures into a single message reporting all errors.
::testing::AssertionResult&
operator&=
(::testing::AssertionResult& old, const ::testing::AssertionResult& val);

// Checks a parameter is marked as defaulted
::testing::AssertionResult
checkDefaultParameterEntry
(const std::string& paramName, const Teuchos::ParameterEntry& param);

// Checks a parameter is marked as used
::testing::AssertionResult
checkUsedParameterEntry
(const std::string& paramName, const Teuchos::ParameterEntry& param);

// Checks a parameter is NOT marked as used
::testing::AssertionResult
checkUnusedParameterEntry
(const std::string& paramName, const Teuchos::ParameterEntry& param);

// Compares whether two parameters are equal. Double are treated specially so
// that NaN values compare equally (to match our use in parameter lists). The
// paramName argument is used in error messages when reporting mismatches.
::testing::AssertionResult
compareParameterEntries
( const std::string& paramName
, const Teuchos::ParameterEntry& param
, const Teuchos::ParameterEntry& validParam
);

// Ensures that the given key is either a list OR missing. If mustExist is true
// missing lists are reported as an error.
::testing::AssertionResult
checkSubList
(const Teuchos::ParameterList& list, const std::string& key, bool mustExist = false);

// Applies an AssertionFunction to every parameter in the ParameterList and
// report the results of combining all these assertions. Intended to be used
// with checkDefaultParameterEntry, checkUsedParameterEntry, and
// checkUnusedParameterEntry
::testing::AssertionResult
checkParameters
( const Teuchos::ParameterList& list
, std::function<::testing::AssertionResult(const std::string&, const Teuchos::ParameterEntry&)> fun
);

// Checks that every parameter in checkList corresponds to a matching entry in
// overrideList. Entries missing from overrideList are compared against the
// defaultList (and checks that these entries in checkList are defaulted).
//
// Entries that exist in defaultList and are missing from checkList are
// reported as missing.
::testing::AssertionResult
checkParameterListAgainstDefaultAndOverrides
( const Teuchos::ParameterList& checkList
, const Teuchos::ParameterList& defaultList
, const Teuchos::ParameterList& overrideList = Teuchos::ParameterList()
);

// Like above, but allows values whose defaults are NaN to change.
::testing::AssertionResult
checkParameterListAgainstDefaultAndOverridesAllowNanChanges
( const Teuchos::ParameterList& checkList
, const Teuchos::ParameterList& defaultList
, const Teuchos::ParameterList& overrideList = Teuchos::ParameterList()
);
#endif
