#ifndef SCOREFUNCTIONS_H
#define SCOREFUNCTIONS_H

#include <functional>

#include "Teuchos_RCP.hpp"

class Epetra_Vector;
class Epetra_MultiVector;

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_default_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3);

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_projected_default_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3,
    Teuchos::RCP<const Epetra_MultiVector> const &V);

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_ocean_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3);

std::function<double(Teuchos::RCP<const Epetra_Vector> const &)>
get_projected_ocean_score_function(
    Teuchos::RCP<const Epetra_Vector> const &sol1,
    Teuchos::RCP<const Epetra_Vector> const &sol2,
    Teuchos::RCP<const Epetra_Vector> const &sol3,
    Teuchos::RCP<const Epetra_MultiVector> const &V);

#endif
