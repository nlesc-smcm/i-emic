#include "TestDefinitions.H"

#include "Combined_MultiVec.H"
#include "TRIOS_Domain.H"
#include "Utils.H"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//------------------------------------------------------------------
namespace
{
    Teuchos::RCP<TRIOS::Domain> domain1, domain2;
    Teuchos::RCP<Epetra_Comm> comm;
    Teuchos::RCP<Epetra_Map> map1, map2, map3;
}

//------------------------------------------------------------------
TEST(Combined_MultiVec, Setup)
{
    bool failed = false;
    try
    {
        // create an arbitrary domain
        domain1 = Teuchos::rcp(new TRIOS::Domain(6, 6, 6, 3, 0, 1, 0,
                                                 1, false, 1.0, 1.0, comm));

        domain1->Decomp2D();
        domain2 = Teuchos::rcp(new TRIOS::Domain(8, 4, 4, 2, 0, 1, 0,
                                                 1, false, 1.0, 1.0, comm));
        domain2->Decomp2D();

        map1 = domain1->GetStandardMap();
        map2 = domain1->GetStandardSurfaceMap();
        map3 = domain2->GetStandardMap();

    }
    catch (...)
    {
        failed = true;
    }

    EXPECT_EQ(failed, false);

}

//------------------------------------------------------------------
TEST(Combined_MultiVec, Construction)
{
    bool failed = false;
    try
    {
        Epetra_MultiVector vec1(*map1, 1);
        Epetra_MultiVector vec2(*map2, 1);
        Epetra_MultiVector vec3(*map3, 1);

        // Constructors using maps
        Combined_MultiVec  two_vec(*map1, *map2, 1);
        Combined_MultiVec  three_vec(*map1, *map2, *map3, 1);

        EXPECT_EQ(two_vec.GlobalLength(),
                  vec1.GlobalLength() + vec2.GlobalLength());

        EXPECT_EQ(three_vec.GlobalLength(),
                  vec1.GlobalLength() + vec2.GlobalLength() + vec3.GlobalLength());

        // Randomize contents
        vec1.Random();
        vec2.Random();
        vec3.Random();

        // Unwanted construction
        Epetra_MultiVector vec4(*map3, 4); // vec4 has a different
                                           // number of Epetra_Vectors
        EXPECT_DEATH(Combined_MultiVec(vec1, vec4), ".*Assertion.*failed.*");
        EXPECT_DEATH(Combined_MultiVec(vec1, vec3, vec4), ".*Assertion.*failed.*");

        // Constructors using multivectors
        two_vec   = Combined_MultiVec(vec1, vec2);
        three_vec = Combined_MultiVec(vec1, vec2, vec3);

        double norm1 = Utils::norm(vec1);
        double norm2 = Utils::norm(vec2);
        double norm3 = Utils::norm(vec3);

        double norm_two_vec   = Utils::norm(two_vec);
        double norm_three_vec = Utils::norm(three_vec);

        std::cout << "||vec1|| = " << norm1 << std::endl;
        std::cout << "||vec2|| = " << norm2 << std::endl;
        std::cout << "||vec3|| = " << norm3 << std::endl;

        EXPECT_EQ(norm_two_vec,
                  sqrt(pow(norm1,2)+pow(norm2,2)));

        EXPECT_EQ(norm_three_vec,
                  sqrt(pow(norm1,2)+pow(norm2,2)+pow(norm3,2)));

        // Deep copy through construction
        Combined_MultiVec copy1(two_vec);
        EXPECT_EQ( Utils::norm(copy1), Utils::norm(two_vec) );

        // Check that it is indeed a copy
        copy1.Scale(2.0);
        EXPECT_NE( Utils::norm(copy1), Utils::norm(two_vec) );

        // Now for three multivectors
        Combined_MultiVec copy2(three_vec);
        EXPECT_EQ( Utils::norm(copy2), Utils::norm(three_vec) );

        // Check that it is indeed a copy
        copy2.Scale(2.0);
        EXPECT_NE( Utils::norm(copy2), Utils::norm(three_vec) );

        // Check assignment
        two_vec   = copy1;
        three_vec = copy2;
        EXPECT_EQ( Utils::norm(copy2), Utils::norm(three_vec) );
        EXPECT_EQ( Utils::norm(copy1), Utils::norm(two_vec)   );

        // Let three Epetra_MultiVectors contain 10 Epetra_Vectors
        Combined_MultiVec three_vec_ten(*map1, *map2, *map3, 10);

        // Randomize contents
        three_vec_ten.SetSeed(99);
        three_vec_ten.Random();
        double originalNorm = Utils::norm(three_vec_ten);

        // Check norms
        std::vector<double> twoNorms_ten(10);
        std::vector<double> infNorms_ten(10);
        std::vector<double> oneNorms_ten(10);

        three_vec_ten.Norm2(  twoNorms_ten);
        three_vec_ten.NormInf(infNorms_ten);
        three_vec_ten.Norm1(  oneNorms_ten);

        for (auto &el: twoNorms_ten)
            EXPECT_NE(el, 0.0);

        for (auto &el: infNorms_ten)
            EXPECT_NE(el, 0.0);

        for (auto &el: oneNorms_ten)
            EXPECT_NE(el, 0.0);

        // Create new Combined_MultiVec with a subset of the above
        std::vector<int> index = {0,2,3,5,9};
        Combined_MultiVec three_vec_five(Copy, three_vec_ten, index);

        std::vector<double> infNorms_five(5);
        three_vec_five.NormInf(infNorms_five);
        EXPECT_EQ(infNorms_five[3], infNorms_ten[5]);

        // Check all norms of three_vec_ten
        double twoNorm, infNorm, oneNorm, tmp;
        for (int v = 0; v != 10; ++v)
        {
            twoNorm = 0.0;
            infNorm = 0.0;
            oneNorm = 0.0;

            for (int i = 0; i != three_vec_ten.Size(); ++i)
            {
                (*three_vec_ten(i))(v)->Norm2(&tmp);
                twoNorm += pow(tmp,2);
                (*three_vec_ten(i))(v)->NormInf(&tmp);
                infNorm = std::max(infNorm, tmp);
                (*three_vec_ten(i))(v)->Norm1(&tmp);
                oneNorm += tmp;
            }

            EXPECT_EQ(oneNorms_ten[v], oneNorm);
            EXPECT_EQ(infNorms_ten[v], infNorm);
            EXPECT_EQ(twoNorms_ten[v], sqrt(twoNorm));
        }

        // Check whether dataAccess = View does give a view of a
        // subset of vectors.
        Combined_MultiVec view(View, three_vec_ten, index);
        std::vector<double> twoNorms_five_view(5);
        std::vector<double> twoNorms_five_copy(5);
        view.Norm2(twoNorms_five_view);           // view
        three_vec_five.Norm2(twoNorms_five_copy); // copy
        EXPECT_EQ((int) twoNorms_five_view.size(), view.NumVectors());
        for (int i = 0; i != view.NumVectors(); ++i)
        {
            EXPECT_EQ(twoNorms_five_view[i], twoNorms_ten[index[i]]);
        }

        // Randomize the original vector again
        three_vec_ten.Random();
        view.Norm2(twoNorms_five_view);
        three_vec_ten.Norm2(twoNorms_ten);
        three_vec_five.Norm2(twoNorms_five_copy);
        for (int i = 0; i != view.NumVectors(); ++i)
        {
            EXPECT_EQ(twoNorms_five_view[i], twoNorms_ten[index[i]]);
            EXPECT_NE(twoNorms_five_view[i],
                      twoNorms_five_copy[i]);
        }

        // Set seed again contents
        three_vec_ten.SetSeed(99);
        three_vec_ten.Random();
        EXPECT_EQ(Utils::norm(three_vec_ten), originalNorm);
    }
    catch (...)
    {
        failed = true;
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

    // Get rid of possibly parallel objects for a clean ending.
    domain1 = Teuchos::null;
    domain2 = Teuchos::null;

    comm->Barrier();
    std::cout << "TEST exit code proc #" << comm->MyPID()
              << " " << out << std::endl;

    MPI_Finalize();
    return out;
}
