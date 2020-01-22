#include "ScoreFunctions.H"

#include "ProjectedThetaModel.H"

#include "Teuchos_RCP.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

double norm2(Teuchos::RCP<const Epetra_Vector> const &vec)
{
    double nrm;
    CHECK_ZERO(vec->Norm2(&nrm));
    return nrm;
}

Teuchos::RCP<Epetra_MultiVector> dot(
    Epetra_MultiVector const &x, Epetra_MultiVector const &y)
{
    TIMER_SCOPE("ProjectedThetaModel: dot");
    int m = x.NumVectors();
    int n = y.NumVectors();

    Epetra_LocalMap map(m, 0, x.Comm());
    Teuchos::RCP<Epetra_MultiVector> out = Teuchos::rcp(new Epetra_MultiVector(map, n));

    CHECK_ZERO(out->Multiply('T', 'N', 1.0, x, y, 0.0));
    return out;
}

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_default_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3)
{
    INFO("Using the default score function");

    Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp(new Epetra_Vector(*sol1));
    CHECK_ZERO(diff->Update(-1.0, *sol2, 1.0));
    double nrm = norm2(diff);

    double dist_factor = 0.5;
    if (sol3 != Teuchos::null)
    {
        Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp(new Epetra_Vector(*sol1));
        CHECK_ZERO(diff->Update(-1.0, *sol3, 1.0));
        dist_factor = norm2(diff) / nrm;
    }
    INFO("distance factor = " << dist_factor);

    return [nrm, dist_factor, sol1, sol2](
        Teuchos::RCP<const Epetra_Vector> const &x) {
        Teuchos::RCP<Epetra_Vector> d1v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d1v->Update(-1.0, *sol1, 1.0));
        double d1 = norm2(d1v) / nrm;
        Teuchos::RCP<Epetra_Vector> d2v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d2v->Update(-1.0, *sol2, 1.0));
        double d2 = norm2(d2v) / nrm;
        double dist = dist_factor - dist_factor * exp(-0.5 * pow(d1 / 0.25, 2.))
            + (1.0 - dist_factor) * exp(-0.5 * pow(d2 / 0.25, 2.));
        INFO("distance = " << dist);
        return dist;
    };
}

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_projected_default_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3,
    Teuchos::RCP<const Epetra_MultiVector> const &V)
{
    INFO("Using the default projected score function");

    auto VV = dot(*V, *V);

    auto projected_norm = [VV](Teuchos::RCP<const Epetra_Vector> const &x) {
        auto tmp = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(tmp->Multiply('N', 'N', 1.0, *VV, *x, 0.0));
        return sqrt((*dot(*x, *tmp))[0][0]);
    };

    auto diff = Teuchos::rcp(new Epetra_Vector(*sol1));
    CHECK_ZERO(diff->Update(-1.0, *sol2, 1.0));
    double nrm = projected_norm(diff);

    double dist_factor = 0.5;
    if (sol3 != Teuchos::null)
    {
        auto diff = Teuchos::rcp(new Epetra_Vector(*sol1));
        CHECK_ZERO(diff->Update(-1.0, *sol3, 1.0));
        double nrm2 = projected_norm(diff);
        dist_factor = nrm2 / nrm;
    }
    INFO("distance factor = " << dist_factor);

    return [nrm, dist_factor, projected_norm, sol1, sol2](
        Teuchos::RCP<const Epetra_Vector> const &x) {
        Teuchos::RCP<Epetra_Vector> d1v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d1v->Update(-1.0, *sol1, 1.0));
        double d1 =  projected_norm(d1v) / nrm;
        Teuchos::RCP<Epetra_Vector> d2v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d2v->Update(-1.0, *sol2, 1.0));
        double d2 =  projected_norm(d2v) / nrm;
        double dist = dist_factor - dist_factor * exp(-0.5 * pow(d1 / 0.25, 2.))
            + (1.0 - dist_factor) * exp(-0.5 * pow(d2 / 0.25, 2.));
        INFO("distance = " << dist);
        return dist;
    };
}

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_ocean_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3)
{
    INFO("Using the ocean score function");

    int vvar = 1;
    int dof = 6;

    auto var_norm = [dof](Teuchos::RCP<const Epetra_Vector> const &x, int var) {
        double tmp = 0;
        int n = x->MyLength();
        auto &map = x->Map();
        Epetra_Vector const &xref = *x;
        for (int i = 0; i < n; i++)
            if (map.GID(i) % dof == var)
                tmp += xref[i] * xref[i];
        double out = 0;
        if (x->DistributedGlobal())
        {
            CHECK_ZERO(x->Comm().SumAll(&tmp, &out, 1));
        }
        else
        {
            out = tmp;
        }
        return sqrt(out);
    };

    std::vector<double> nrm(dof);
    for (int var = 0; var < dof; var++)
    {
        Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp(new Epetra_Vector(*sol1));
        CHECK_ZERO(diff->Update(-1.0, *sol2, 1.0));
        nrm[var] = var_norm(diff, var);
    }

    double dist_factor = 0.5;
    if (sol3 != Teuchos::null)
    {
        Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp(new Epetra_Vector(*sol1));
        CHECK_ZERO(diff->Update(-1.0, *sol3, 1.0));
        double nrm2 = var_norm(diff, vvar);
        dist_factor = nrm2 / nrm[vvar];
    }
    INFO("distance factor = " << dist_factor);

    return [nrm, dist_factor, var_norm, vvar, sol1, sol2](
        Teuchos::RCP<const Epetra_Vector> const &x) {
        // // debug
        // for (int var = 0; var < dof; var++)
        // {
        //     Teuchos::RCP<const Epetra_Vector> d1v = Teuchos::rcp(new Epetra_Vector(*x));
        //     CHECK_ZERO(d1v->Update(-1.0, *sol1, 1.0));
        //     double d1 =  var_norm(d1v, var) / nrm[var];
        //     Teuchos::RCP<const Epetra_Vector> d2v = Teuchos::rcp(new Epetra_Vector(*x));
        //     CHECK_ZERO(d2v->Update(-1.0, *sol2, 1.0));
        //     double d2 =  var_norm(d2v, var) / nrm[var];
        //     double dist = dist_factor - dist_factor * exp(-0.5 * pow(d1 / 0.25, 2.))
        //     + (1.0 - dist_factor) * exp(-0.5 * pow(d2 / 0.25, 2.));
        //     INFO("distance " << var << " = " << dist << ", " << d1  << ", " << d2);
        // }
        // // end debug
        Teuchos::RCP<Epetra_Vector> d1v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d1v->Update(-1.0, *sol1, 1.0));
        double d1 =  var_norm(d1v, vvar) / nrm[vvar];
        Teuchos::RCP<Epetra_Vector> d2v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d2v->Update(-1.0, *sol2, 1.0));
        double d2 =  var_norm(d2v, vvar) / nrm[vvar];
        double dist = dist_factor - dist_factor * exp(-0.5 * pow(d1 / 0.25, 2.))
            + (1.0 - dist_factor) * exp(-0.5 * pow(d2 / 0.25, 2.));
        INFO("distance = " << dist);
        return dist;
    };
}

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_projected_ocean_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3,
    Teuchos::RCP<const Epetra_MultiVector> const &V)
{
    INFO("Using the projected ocean score function");

    int vvar = 1;
    int dof = 6;

    auto vvec = Teuchos::rcp(new Epetra_Vector(V->Map()));
    CHECK_ZERO(vvec->PutScalar(0.0));
    for (int i = 0; i < vvec->MyLength(); i++)
        if (vvec->Map().GID(i) % dof == vvar)
            (*vvec)[i] = 1.0;

    auto Vvvec = Teuchos::rcp(new Epetra_MultiVector(*V));
    CHECK_ZERO(Vvvec->Multiply(1.0, *vvec, *V, 0.0));

    auto VvvV = dot(*Vvvec, *Vvvec);

    auto projected_v_norm = [VvvV](Teuchos::RCP<const Epetra_Vector> const &x) {
        auto tmp = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(tmp->Multiply('N', 'N', 1.0, *VvvV, *x, 0.0));
        return sqrt((*dot(*x, *tmp))[0][0]);
    };

    auto diff = Teuchos::rcp(new Epetra_Vector(*sol1));
    CHECK_ZERO(diff->Update(-1.0, *sol2, 1.0));
    double nrm = projected_v_norm(diff);

    double dist_factor = 0.5;
    if (sol3 != Teuchos::null)
    {
        auto diff = Teuchos::rcp(new Epetra_Vector(*sol1));
        CHECK_ZERO(diff->Update(-1.0, *sol3, 1.0));
        double nrm2 = projected_v_norm(diff);
        dist_factor = nrm2 / nrm;
    }
    INFO("distance factor = " << dist_factor);

    return [nrm, dist_factor, projected_v_norm, sol1, sol2](
        Teuchos::RCP<const Epetra_Vector> const &x) {
        Teuchos::RCP<Epetra_Vector> d1v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d1v->Update(-1.0, *sol1, 1.0));
        double d1 =  projected_v_norm(d1v) / nrm;
        Teuchos::RCP<Epetra_Vector> d2v = Teuchos::rcp(new Epetra_Vector(*x));
        CHECK_ZERO(d2v->Update(-1.0, *sol2, 1.0));
        double d2 =  projected_v_norm(d2v) / nrm;
        double dist = dist_factor - dist_factor * exp(-0.5 * pow(d1 / 0.25, 2.))
        + (1.0 - dist_factor) * exp(-0.5 * pow(d2 / 0.25, 2.));
        INFO("distance = " << dist);
        return dist;
    };
}
