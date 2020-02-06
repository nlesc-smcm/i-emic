#include "TestDefinitions.H"

//------------------------------------------------------------------
TEST(ParameterList, UnusedParameters)
{
    Teuchos::ParameterList list("Test List");

    list.set("Test", true);
    EXPECT_TRUE(checkParameters(list, checkUnusedParameterEntry));
    EXPECT_FALSE(checkParameters(list, checkUsedParameterEntry));
}

//------------------------------------------------------------------
TEST(ParameterList, UnusedSublistParameters)
{
    Teuchos::ParameterList list("Test List");

    list.sublist("Sublist").set("Test", true);
    EXPECT_TRUE(checkParameters(list, checkUnusedParameterEntry));
    EXPECT_FALSE(checkParameters(list, checkUsedParameterEntry));
}

//------------------------------------------------------------------
TEST(ParameterList, DefaultParameters)
{
    Teuchos::ParameterList list("Test List");

    list.get("Test", true);
    EXPECT_TRUE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
TEST(ParameterList, NotDefaultParameters)
{
    Teuchos::ParameterList list("Test List");

    list.set("Test", true);
    EXPECT_FALSE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
TEST(ParameterList, DefaultSublistParameters)
{
    Teuchos::ParameterList list("Test List");

    list.sublist("Sublist").get("Test", true);
    EXPECT_TRUE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
TEST(ParameterList, NotDefaultSublistParameters)
{
    Teuchos::ParameterList list("Test List");

    list.sublist("Sublist").set("Test", true);
    EXPECT_FALSE(checkParameters(list, checkDefaultParameterEntry));
}

//------------------------------------------------------------------
TEST(ParameterList, DefaultEqualsSelf)
{
    Teuchos::ParameterList list("Test List");

    list.get("Test", true);
    list.sublist("Sublist").get("Test", true);
    EXPECT_TRUE(checkParameterList(list, list));
}

//------------------------------------------------------------------
TEST(ParameterList, DefaultEqualsSelf2)
{
    Teuchos::ParameterList list("Test List");

    list.get("Test", true);
    list.sublist("Sublist").get("Test", true);
    EXPECT_TRUE(checkParameterList(list, list, list));
}

//------------------------------------------------------------------
TEST(ParameterList, CopiesMatch)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);

    EXPECT_TRUE(checkParameterList(list, defaultList));
}


//------------------------------------------------------------------
TEST(ParameterList, NoMissingParameters)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);
    list.remove("Test");

    EXPECT_FALSE(checkParameterList(list, defaultList));
}

//------------------------------------------------------------------
TEST(ParameterList, NoMissingSublist)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);
    list.remove("Sublist");

    EXPECT_FALSE(checkParameterList(list, defaultList));
}

//------------------------------------------------------------------
TEST(ParameterList, MissingOverriddenParameter)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList overrideList("Override List");
    overrideList.set("Test", false);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);

    EXPECT_FALSE(checkParameterList(list, defaultList, overrideList));
}

//------------------------------------------------------------------
TEST(ParameterList, MissingOverriddenSublist)
{
    Teuchos::ParameterList defaultList("Default List");

    defaultList.get("Test", true);
    defaultList.sublist("Sublist").get("Test", true);

    Teuchos::ParameterList overrideList("Override List");
    overrideList.sublist("Sublist").set("Test", false);

    Teuchos::ParameterList list("Test List");
    list.setParameters(defaultList);

    EXPECT_FALSE(checkParameterList(list, defaultList, overrideList));
}

//------------------------------------------------------------------
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

    EXPECT_TRUE(checkParameterList(list, defaultList, overrideList));
}

//------------------------------------------------------------------
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

    EXPECT_TRUE(checkParameterList(list, defaultList, overrideList));
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
