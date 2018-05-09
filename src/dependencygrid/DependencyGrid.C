#include "DependencyGrid.H"

//=============================================================================
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//                                                                           //
//                                                                           //
// DependencyGrid implementation                                             //
//                                                                           //
//                                                                           //
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//=============================================================================
DependencyGrid::DependencyGrid(int n, int m, int l, int np, int nun)
    :
    grid_(n, m, l, np, nun, nun),
    n_(n),
    m_(m),
    l_(l),
    np_(np),
    nun_(nun)
{}

//-----------------------------------------------------------------------------
DependencyGrid::~DependencyGrid()
{}

//-----------------------------------------------------------------------------
double &DependencyGrid::operator() (int i, int j, int k, int loc, int A, int B)
{
    // converting to 0-based
    return grid_(i-1, j-1, k-1, loc-1, A-1, B-1);
}

//-----------------------------------------------------------------------------
double DependencyGrid::get(int i, int j, int k, int loc, int A, int B)
{
    // converting to 0-based
    return grid_(i-1, j-1, k-1, loc-1, A-1, B-1);
}

//-----------------------------------------------------------------------------
void DependencyGrid::set(int i, int j, int k, int loc, int A, int B, double value)
{
    // converting to 0-based
    grid_(i-1, j-1, k-1, loc-1, A-1, B-1) = value;
}

//-----------------------------------------------------------------------------
void DependencyGrid::set(int const (&range)[8], int A, int B, double value)
{
    for (int i = range[0]; i != range[1]+1; ++i)
        for (int j = range[2]; j != range[3]+1; ++j)
            for (int k = range[4]; k != range[5]+1; ++k)
                for (int loc = range[6]; loc != range[7]+1; ++loc)
                {
                    // converting to 0-based
                    grid_(i-1, j-1, k-1, loc-1, A-1, B-1) = value;
                }
}

//-----------------------------------------------------------------------------
void DependencyGrid::set(int const (&range)[8], int A, int B, Atom &atom)
{
    for (int i = range[0]; i != range[1]+1; ++i)
        for (int j = range[2]; j != range[3]+1; ++j)
            for (int k = range[4]; k != range[5]+1; ++k)
                for (int loc = range[6]; loc != range[7]+1; ++loc)
                {
                    // converting to 0-based
                    grid_(i-1, j-1, k-1, loc-1, A-1, B-1) =
                        atom.get(i,j,k,loc);
                }
}

//-----------------------------------------------------------------------------
void DependencyGrid::zero()
{
    grid_.assign(0.0);
}


//=============================================================================
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//                                                                           //
//                                                                           //
// Atom implementation                                                       //
//                                                                           //
//                                                                           //
// / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / //
//=============================================================================
Atom::Atom(int n, int m, int l, int np)
    :
    atom_(n, m, l, np),
    n_(n),
    m_(m),
    l_(l),
    np_(np)
{}

//-----------------------------------------------------------------------------
Atom::~Atom()
{}

//-----------------------------------------------------------------------------
double Atom::get(int i, int j, int k, int loc)
{
    // converting to 0-based
    return atom_(i-1, j-1, k-1, loc-1);
}

//-----------------------------------------------------------------------------
void Atom::set(int i, int j, int k, int loc, double value)
{
    // converting to 0-based
    atom_(i-1, j-1, k-1, loc-1) = value;
}

//-----------------------------------------------------------------------------
void Atom::set(int const (&range)[6], int loc, double value)
{
    for (int i = range[0]; i != range[1]+1; ++i)
        for (int j = range[2]; j != range[3]+1; ++j)
            for (int k = range[4]; k != range[5]+1; ++k)
            {
                atom_(i-1, j-1, k-1, loc-1) = value;
            }
}

//-----------------------------------------------------------------------------
// this = scalThis*this+scalA*A+scalB*B
void Atom::update(double scalarThis,
                  double scalarA, Atom &A,
                  double scalarB, Atom &B)
{
    for (int i = 1; i != n_+1; ++i)
        for (int j = 1; j != m_+1; ++j)
            for (int k = 1; k != l_+1; ++k)
                for (int loc = 1; loc != np_+1; ++loc)
                {
                    // converting to 0-based
                    atom_(i-1, j-1, k-1, loc-1) =
                        scalarThis * atom_(i-1, j-1, k-1, loc-1) +
                        scalarA*A.get(i, j, k, loc) +
                        scalarB*B.get(i, j, k, loc);
                }
}

//-----------------------------------------------------------------------------
// this = scalThis*this+scalA*A+scalB*B+scalC*C
void Atom::update(double scalarThis,
                  double scalarA, Atom &A,
                  double scalarB, Atom &B,
                  double scalarC, Atom &C)
{
    for (int i = 1; i != n_+1; ++i)
        for (int j = 1; j != m_+1; ++j)
            for (int k = 1; k != l_+1; ++k)
                for (int loc = 1; loc != np_+1; ++loc)
                {
                    // converting to 0-based
                    atom_(i-1, j-1, k-1, loc-1) =
                        scalarThis * atom_(i-1, j-1, k-1, loc-1) +
                        scalarA*A.get(i, j, k, loc) +
                        scalarB*B.get(i, j, k, loc) +
                        scalarC*C.get(i, j, k, loc);
                }
}

//-----------------------------------------------------------------------------
// this = scalThis*this
void Atom::scale(double scalarThis)
{
    for (int i = 1; i != n_+1; ++i)
        for (int j = 1; j != m_+1; ++j)
            for (int k = 1; k != l_+1; ++k)
                for (int loc = 1; loc != np_+1; ++loc)
                {
                    // converting to 0-based
                    atom_(i-1, j-1, k-1, loc-1) =
                        scalarThis * atom_(i-1, j-1, k-1, loc-1);
                }
}
