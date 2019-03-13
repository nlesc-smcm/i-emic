//-----------------------------------------------------------------------------
/*
   Dependency grid:
    A multidimensional array describing the dependencies among the unknowns:
    grid(i,j,k,21,U,V) = c   <=>   (...) * d/dt U|(i,j,k) = ... + c * V|(i-1,j+1,k+1)
    grid(i,j,k,13,U,V) = d   <=>   (...) * d/dt U|(i,j,k) = ... + d * V|(i,j-1,k-1)

*/
//-----------------------------------------------------------------------------

#ifndef DEPENDENCYGRID_H
#define DEPENDENCYGRID_H

#include "MultiArray.H"

class Atom;

class DependencyGrid
{
    MultiArray<double, 6> grid_;
    int n_, m_, l_, np_, nun_;

public:
    DependencyGrid(int n, int m, int l, int np, int nun);
    ~DependencyGrid();

    double &operator() (int i, int j, int k, int loc, int A, int B);

    // these get and set members are all 1-based!
    double get(int i, int j, int k, int loc, int A, int B);
    void   set(int i, int j, int k, int loc, int A, int B, double value);
    void   set(int const (&range)[8], int A, int B, Atom &atom);
    void   set(int const (&range)[8], int A, int B, double value);
    void   set(int i, int j, int k, int loc, int const (&range)[4], double value);

    void   add(double scalar, Atom &atom);
    void   zero();

};

//-----------------------------------------------------------------------------
/*
//! Atom:
//! A multidimensional array describing anonymous dependencies among neighbours:
//! atom(i,j,k,21) = c   <=>   (...) * d/dt {}|(i,j,k) = ... + c * {}|(i-1,j+1,k+1)
//! atom(i,j,k,13) = d   <=>   (...) * d/dt {}|(i,j,k) = ... + d * {}|(i,j-1,k-1)
*/
//-----------------------------------------------------------------------------
class Atom
{
    MultiArray<double, 4> atom_;
    int n_, m_, l_, np_;

public:
    Atom(int n, int m, int l, int np);
    ~Atom();

    double get(int i, int j, int k, int loc);
    void   set(int i, int j, int k, int loc, double value);
    void   set(int const (&range)[6], int loc, double value);

    void update(double scalarThis,
                double scalarA, Atom &A,
                double scalarB, Atom &B);

    void update(double scalarThis,
                double scalarA, Atom &A,
                double scalarB, Atom &B,
                double scalarC, Atom &C);

    void scale(double scalarThis);

    // this = vec.*this (pointwise) along dimension dim
    void multiply(int dim, std::vector<double> &vec, double scalarThis);

};


#endif
