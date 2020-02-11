#include "TestDefinitions.H"

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

//------------------------------------------------------------------
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // -------------------------------------------------------
    // TESTING
    int out = RUN_ALL_TESTS();
    // -------------------------------------------------------

    return out;
}
