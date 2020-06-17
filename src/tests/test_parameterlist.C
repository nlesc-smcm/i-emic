#include "TestDefinitions.H"

#include "Continuation.H"
#include "Ocean.H"
#include "THCM.H"

namespace // local unnamed namespace (similar to static in C)
{
    Teuchos::RCP<Epetra_Comm>  comm;
}

//------------------------------------------------------------------
// Check that parameters that are explicitly set are reported as unused
TEST(ParameterList, UnusedParameters)
{
    Teuchos::ParameterList list("Test List");

    list.set("Test", true);
    EXPECT_TRUE(checkParameters(list, checkUnusedParameterEntry));
    EXPECT_FALSE(checkParameters(list, checkUsedParameterEntry));
}

//------------------------------------------------------------------
// Check that parameters that are explicitly set in sublists are reported as
// unused
TEST(ParameterList, UnusedSublistParameters)
{
    Teuchos::ParameterList list("Test List");

    list.sublist("Sublist").set("Test", true);
    EXPECT_TRUE(checkParameters(list, checkUnusedParameterEntry));
    EXPECT_FALSE(checkParameters(list, checkUsedParameterEntry));
}

//------------------------------------------------------------------
// Check that parameters set via get defaulting are marked as defaulted
TEST(ParameterList, DefaultParameters)
{
    Teuchos::ParameterList list("Test List");

    list.get("Test", true);
    EXPECT_TRUE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
// Check that set parameters are NOT marked as defaulted
TEST(ParameterList, NotDefaultParameters)
{
    Teuchos::ParameterList list("Test List");

    list.set("Test", true);
    EXPECT_FALSE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
// Check that sublist parameters set via get defaulting are marked as defaulted
TEST(ParameterList, DefaultSublistParameters)
{
    Teuchos::ParameterList list("Test List");

    list.sublist("Sublist").get("Test", true);
    EXPECT_TRUE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
// Check that set sublist parameters set are NOT marked as defaulted
TEST(ParameterList, NotDefaultSublistParameters)
{
    Teuchos::ParameterList list("Test List");

    list.sublist("Sublist").set("Test", true);
    EXPECT_FALSE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
// Check that a default list matches itself
TEST(ParameterList, DefaultEqualsSelf)
{
    Teuchos::ParameterList list("Test List");

    list.get("Test", true);
    list.sublist("Sublist").get("Test", true);
    // checkParameterListAgainstDefaultAndOverrides(list1, list2) checks that
    // every element in list1 matches the elements in list2 and is marked as
    // defaulted
    EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(list, list));
}

//------------------------------------------------------------------
// Check that a default list with sublists matches itself
TEST(ParameterList, DefaultEqualsSelf2)
{
    Teuchos::ParameterList list("Test List");

    list.get("Test", true);
    list.sublist("Sublist").get("Test", true);
    // checkParameterListAgainstDefaultAndOverrides(list1, list2, list3) checks
    // that every element in list1 matches the elements in list3. If list3 is
    // missing an element it is matched with list2 and must be marked as
    // defaulted
    EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(list, list, list));
}

//------------------------------------------------------------------
// Check that the copy of a list matches the original default list
TEST(ParameterList, CopiesMatch)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);

    // checkParameterListAgainstDefaultAndOverrides(list1, list2) checks that
    // every element in list1 matches the elements in list2 and is marked as
    // defaulted
    EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(list, defaultList));
}

//------------------------------------------------------------------
// Check that removing a parameter reports an error
TEST(ParameterList, NoMissingParameters)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);
    list.remove("Test");

    // checkParameterListAgainstDefaultAndOverrides(list1, list2) checks that
    // every element in list1 matches the elements in list2 and is marked as
    // defaulted
    EXPECT_FALSE(checkParameterListAgainstDefaultAndOverrides(list, defaultList));
}

//------------------------------------------------------------------
// Check that removing a sublist reports an error
TEST(ParameterList, NoMissingSublist)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);
    list.remove("Sublist");

    // checkParameterListAgainstDefaultAndOverrides(list1, list2) checks that
    // every element in list1 matches the elements in list2 and is marked as
    // defaulted
    EXPECT_FALSE(checkParameterListAgainstDefaultAndOverrides(list, defaultList));
}

//------------------------------------------------------------------
// Check that parameters in the override list must be present
TEST(ParameterList, MissingOverriddenParameter)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList overrideList("Override List");
    overrideList.set("Test", false);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);

    // checkParameterListAgainstDefaultAndOverrides(list1, list2, list3) checks
    // that every element in list1 matches the elements in list3. If list3 is
    // missing an element it is matched with list2 and must be marked as
    // defaulted
    EXPECT_FALSE(checkParameterListAgainstDefaultAndOverrides(list, defaultList, overrideList));
}

//------------------------------------------------------------------
// Check that parameters in sublists in the override list must be present
TEST(ParameterList, MissingOverriddenSublist)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList overrideList("Override List");
    overrideList.sublist("Sublist").set("Test", false);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);

    // checkParameterListAgainstDefaultAndOverrides(list1, list2, list3) checks
    // that every element in list1 matches the elements in list3. If list3 is
    // missing an element it is matched with list2 and must be marked as
    // defaulted
    EXPECT_FALSE(checkParameterListAgainstDefaultAndOverrides(list, defaultList, overrideList));
}

//------------------------------------------------------------------
// Check that a combination of overriden parameters and defaulted parameters
// validates
TEST(ParameterList, MatchOverriddenParameter)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList overrideList("Override List");
    overrideList.set("Test", false);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);
    list.set("Test", false);

    // checkParameterListAgainstDefaultAndOverrides(list1, list2, list3) checks
    // that every element in list1 matches the elements in list3. If list3 is
    // missing an element it is matched with list2 and must be marked as
    // defaulted
    EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(list, defaultList, overrideList));
}

//------------------------------------------------------------------
// Check that a combination of overriden parameters in sublist and defaulted
// parameters validates
TEST(ParameterList, MatchOverriddenSublist)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList overrideList("Override List");
    overrideList.sublist("Sublist").set("Test", false);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);
    list.sublist("Sublist").set("Test", false);

    // checkParameterListAgainstDefaultAndOverrides(list1, list2, list3) checks
    // that every element in list1 matches the elements in list3. If list3 is
    // missing an element it is matched with list2 and must be marked as
    // defaulted
    EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(list, defaultList, overrideList));
}

TEST(ParameterList, CheckNanChangeFailure)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", std::numeric_limits<double>::quiet_NaN());

    Teuchos::ParameterList list("Test List");
    list.set("Test", 15.0);

    // checkParameterListAgainstDefaultAndOverrides(list1, list2, list3) checks
    // that every element in list1 matches the elements in list3. If list3 is
    // missing an element it is matched with list2 and must be marked as
    // defaulted
    EXPECT_FALSE(checkParameterListAgainstDefaultAndOverrides(list, defaultList));
}

TEST(ParameterList, CheckNanChangeSuccess)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", std::numeric_limits<double>::quiet_NaN());

    Teuchos::ParameterList list("Test List");
    list.get("Test", 15.0);

    // checkParameterListAgainstDefaultAndOverrides(list1, list2, list3) checks
    // that every element in list1 matches the elements in list3. If list3 is
    // missing an element it is matched with list2 and must be marked as
    // defaulted
    EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(list, defaultList));
}

//------------------------------------------------------------------
TEST(THCMParameterList, DefaultInitialization)
{
    Teuchos::RCP<THCM> thcm;

    bool failed = false;
    try
    {
        // Copy of the default parameters
        const Teuchos::ParameterList defaultParams = THCM::getDefaultInitParameters();

        // Empty parameter configuration
        Teuchos::ParameterList thcmParams("Test List");

        ::testing::internal::CaptureStdout();
        thcm = Teuchos::rcp(new THCM(thcmParams, comm));
        ::testing::internal::GetCapturedStdout();

        // Copy of the configuration reported by Ocean
        const Teuchos::ParameterList& currentParams = thcm->getParameters();

        // Don't check parameters are used, as default configuration has some
        // unused variables (the mask files)

        // Check that every parameter in oceanParams matches defaultParams, and
        // check that they're all reported as defaulted
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(thcmParams, defaultParams));
        EXPECT_TRUE(checkParameters(thcmParams, checkDefaultParameterEntry));

        // Check that every parameter reported by Ocean matches defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(currentParams, defaultParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(THCMParameterList, Initialization)
{
    Teuchos::RCP<THCM> thcm;

    thcm = Teuchos::null;
    bool failed = false;
    try
    {
        // Copy of the input parameters to validate against
        const Teuchos::ParameterList startParams(Utils::obtainParams("ocean_params.xml", "THCM parameters")->sublist("THCM"));

        // Copy of the default parameters
        const Teuchos::ParameterList defaultParams(THCM::getDefaultInitParameters());
        // The input parameters to validate with
        Teuchos::ParameterList thcmParams(Utils::obtainParams("ocean_params.xml", "THCM parameters")->sublist("THCM"));

        ::testing::internal::CaptureStdout();
        thcm = Teuchos::rcp(new THCM(thcmParams, comm));
        ::testing::internal::GetCapturedStdout();

        // Parameters currently reported by Ocean
        const Teuchos::ParameterList& currentParams = thcm->getParameters();

        // Check that every parameter is used
        //
        // This MUST be before checkParameterListAgainstDefaultAndOverrides
        // because it marks everything as used...
        EXPECT_TRUE(checkParameters(thcmParams, checkUsedParameterEntry));

        // Check that every entry in oceanParams corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(thcmParams, defaultParams, startParams));

        // Check that every entry reported by Ocean corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(currentParams, defaultParams, startParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(OceanParameterList, DefaultInitialization)
{
    Teuchos::RCP<Ocean> ocean;

    bool failed = false;
    try
    {
        // Copy of the default parameters
        const Teuchos::ParameterList defaultParams = Ocean::getDefaultInitParameters();

        // Empty parameter configuration
        Teuchos::ParameterList oceanParams("Test List");

        ::testing::internal::CaptureStdout();
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
        ::testing::internal::GetCapturedStdout();

        // Copy of the configuration reported by Ocean
        const Teuchos::ParameterList& currentParams = ocean->getParameters();

        // Don't check parameters are used, as default configuration has some
        // unused variables (the mask files)

        // Check that every parameter in oceanParams matches defaultParams, and
        // check that they're all reported as defaulted
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(oceanParams, defaultParams));
        EXPECT_TRUE(checkParameters(oceanParams, checkDefaultParameterEntry));

        // Check that every parameter reported by Ocean matches defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(currentParams, defaultParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(OceanParameterList, Initialization)
{
    Teuchos::RCP<Ocean> ocean;

    bool failed = false;
    try
    {
        // Copy of the input parameters to validate against
        const Teuchos::ParameterList startParams(*Utils::obtainParams("ocean_params.xml", "Ocean parameters"));

        // Copy of the default parameters
        const Teuchos::ParameterList defaultParams(Ocean::getDefaultInitParameters());
        // The input parameters to validate with
        Teuchos::ParameterList oceanParams(*Utils::obtainParams("ocean_params.xml", "Ocean parameters"));

        ::testing::internal::CaptureStdout();
        ocean = Teuchos::rcp(new Ocean(comm, oceanParams));
        ::testing::internal::GetCapturedStdout();

        // Parameters currently reported by Ocean
        const Teuchos::ParameterList& currentParams = ocean->getParameters();

        // Check that every parameter is used
        //
        // This MUST be before checkParameterListAgainstDefaultAndOverrides
        // because it marks everything as used...
        EXPECT_TRUE(checkParameters(oceanParams, checkUsedParameterEntry));

        // Check that every entry in oceanParams corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(oceanParams, defaultParams, startParams));

        // Check that every entry reported by Ocean corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverridesAllowNanChanges(currentParams, defaultParams, startParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
TEST(ContinuationParameterList, Initialization)
{
    ::testing::internal::CaptureStdout();
    Teuchos::RCP<Ocean> ocean = Teuchos::rcp(new Ocean(comm));
    ::testing::internal::GetCapturedStdout();
    Teuchos::RCP<Continuation<Teuchos::RCP<Ocean>>> continuation;

    bool failed = false;
    try
    {
        // Copy of the input parameters to validate against
        Teuchos::ParameterList startParams;

        std::stringstream destID;
        for (int i = 0; i != 10; ++i)
        {
            destID << "destination " << i;
            startParams.set(destID.str(), 1.0);
            destID.str("");
            destID.clear();
        }

        // Copy of the default parameters
        const Teuchos::ParameterList defaultParams = Continuation<Teuchos::RCP<Ocean>>::getDefaultInitParameters();

        // Empty parameter configuration
        Teuchos::ParameterList continuationParams = startParams;

        ::testing::internal::CaptureStdout();
        continuation = Teuchos::rcp(new Continuation<Teuchos::RCP<Ocean>>(ocean, continuationParams));
        ::testing::internal::GetCapturedStdout();

        // Copy of the configuration reported by Continuation
        const Teuchos::ParameterList& currentParams = continuation->getParameters();

        // Check that every entry in continuationParams corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(continuationParams, defaultParams, startParams));

        // Check that every entry reported by Continuation corresponds to the value in
        // startParams, missing entries are compared against defaultParams
        EXPECT_TRUE(checkParameterListAgainstDefaultAndOverrides(currentParams, defaultParams, startParams));
    }
    catch (...)
    {
        failed = true;
        throw;
    }

    EXPECT_EQ(failed, false);
}

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Initialize the environment:
    comm = initializeEnvironment(argc, argv);
    if (outFile == Teuchos::null)
        throw std::runtime_error("ERROR: Specify output streams");

    ::testing::InitGoogleTest(&argc, argv);

    // -------------------------------------------------------
    // TESTING
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    if (comm->MyPID() == 0)
        printProfile();

    MPI_Finalize();
    return out;
}
