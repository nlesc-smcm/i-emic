#include "TestDefinitions.H"

//------------------------------------------------------------------
double norm(Teuchos::RCP<Epetra_Vector> vec)
{
    double nrm;
    vec->Norm2(&nrm);
    return nrm;
}

//------------------------------------------------------------------
double norm(std::shared_ptr<std::vector<double> > vec)
{
    int dim = (int) vec->size();
    int incX = 1;
    int incY = 1;
    double dot;
    dot = ddot_(&dim, &(*vec)[0], &incX, &(*vec)[0], &incY);
    return sqrt(dot);
}

//------------------------------------------------------------------
std::shared_ptr<std::vector<double> >
getGatheredVector(Teuchos::RCP<Epetra_Vector> vec)
{
    // get size
    int dim = (int) vec->GlobalLength();

    // Create gather scheme: map and importer
    Teuchos::RCP<Epetra_BlockMap> gmap =
        Utils::AllGather(vec->Map());
    Teuchos::RCP<Epetra_Import> imp =
        Teuchos::rcp(new Epetra_Import(*gmap, vec->Map()));

    // Create vector to hold the gathered state
    Teuchos::RCP<Epetra_Vector> gathered =
        Teuchos::rcp(new Epetra_Vector(*gmap));

    // Gather the state into <gathered>
    gathered->Import(*vec, *imp, Insert);

    // Create full state for serial atmosphere
    std::shared_ptr< std::vector<double> > gvec =
        std::make_shared<std::vector<double> >(dim, 0.0);

    // Extract gvec from gathered
    gathered->ExtractCopy(&(*gvec)[0], dim);

    return gvec;
}

//------------------------------------------------------------------
// Functions for checking Teuchos::ParameterList's
::testing::AssertionResult&
operator&=
(::testing::AssertionResult& old, const ::testing::AssertionResult& val)
{
    if (old == ::testing::AssertionSuccess()) {
        old = val;
    } else if (val != ::testing::AssertionSuccess()) {
        old = ::testing::AssertionFailure() << old.message() << val.message();
    }
    return old;
}

::testing::AssertionResult
checkDefaultParameterEntry
(const std::string& paramName, const Teuchos::ParameterEntry& param)
{
    if (!param.isDefault()) {
        return ::testing::AssertionFailure() << std::endl
            << paramName << " expected defaulted parameter" << std::endl;
    }

    return ::testing::AssertionSuccess();
}

::testing::AssertionResult
checkUsedParameterEntry
(const std::string& paramName, const Teuchos::ParameterEntry& param)
{
    if (!param.isUsed()) {
        return ::testing::AssertionFailure() << std::endl
            << paramName << " expected used parameter" << std::endl;
    }

    return ::testing::AssertionSuccess();
}

::testing::AssertionResult
checkUnusedParameterEntry
(const std::string& paramName, const Teuchos::ParameterEntry& param)
{
    if (param.isUsed()) {
        return ::testing::AssertionFailure() << std::endl
            << paramName << " expected unused parameter" << std::endl;
    }

    return ::testing::AssertionSuccess();
}

::testing::AssertionResult
compareParameterEntries
( const std::string& paramName
, const Teuchos::ParameterEntry& param
, const Teuchos::ParameterEntry& validParam
)
{
    if (param.isType<double>() && validParam.isType<double>()) {
        if (std::isnan(param.getValue<double>(nullptr)) && std::isnan(validParam.getValue<double>(nullptr))) {
            return ::testing::AssertionSuccess() << "NaNs treated as equal";
        }
    }

    if (param == validParam) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() << std::endl
            << paramName << ":" << std::endl
            << "    Found:    " << param << std::endl
            << "    Expected: " << validParam << std::endl;
    }
}

::testing::AssertionResult
checkSubList
(const Teuchos::ParameterList& list, const std::string& key, bool mustExist)
{
    if (mustExist && !list.isParameter(key)) {
        return ::testing::AssertionFailure() << key << " sublist is missing and must exist.";
    } else if (list.isParameter(key) && !list.isSublist(key)) {
        return ::testing::AssertionFailure() << key << " is a parameter, but not a sublist";
    }

    return ::testing::AssertionSuccess();
}

::testing::AssertionResult
checkParameters
( const Teuchos::ParameterList& list
, std::function<::testing::AssertionResult(const std::string&, const Teuchos::ParameterEntry&)> fun
)
{
    auto result = ::testing::AssertionSuccess();

    for (const auto& pair : list) {
        const std::string& key = pair.first;
        const Teuchos::ParameterEntry& value = pair.second;

        if (value.isList()) {
            result &= checkParameters(list.sublist(key), fun);
        } else {
            result &= fun(list.name() + "->" + key, value);
        }
    }

    return result;
}

::testing::AssertionResult
checkParameterListAgainstDefaultAndOverrides
( const Teuchos::ParameterList& checkList
, const Teuchos::ParameterList& defaultList
, const Teuchos::ParameterList& overrideList
)
{
    auto result = ::testing::AssertionSuccess();
    std::set<std::string> parameters;
    for (const auto& pair : defaultList) {
        parameters.insert(pair.first);
    }

    for (const auto& pair : checkList) {
        const std::string& key = pair.first;
        const Teuchos::ParameterEntry& value = pair.second;

        parameters.erase(key);

        if (value.isList()) {
            result &= checkSubList(overrideList, key);
            result &= checkSubList(defaultList, key, true);

            Teuchos::ParameterList sublist;
            if (overrideList.isSublist(key)) {
                sublist = overrideList.sublist(key);
            }

            result &= checkParameterListAgainstDefaultAndOverrides(checkList.sublist(key), defaultList.sublist(key), sublist);
        } else {
            std::string name = checkList.name() + "->" + key;

            Teuchos::ParameterEntry validVal;
            if (overrideList.isParameter(key)) {
                validVal = overrideList.getEntry(key);
            } else {
                validVal = defaultList.getEntry(key);
                if (!value.isDefault()) {
                    result &= ::testing::AssertionFailure() << "expected " << name << " to be defaulted";
                }
            }

            result &= compareParameterEntries(name, value, validVal);
        }
    }

    if (parameters.size() != 0) {
        auto error = ::testing::AssertionFailure() << std::endl;
        for (auto& param : parameters) {
            error << "Missing parameter: " << checkList.name() + "->" + param;
        }

        result &= error;
    }

    return result;
}
