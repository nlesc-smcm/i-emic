//// ---------------------------------------------------
//// THIS CODE IS BROKEN BUT IT MAY SERVE AS AN EXAMPLE
//// ---------------------------------------------------

#ifndef TEST_H
#define TEST_H

#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>

#include <Epetra_Vector.h>
#include <EpetraExt_RowMatrixOut.h>

#include <Teuchos_RCP.hpp>

#include "THCM.H"

class SolverInterface : public NOX::Epetra::Interface::Required,
						public NOX::Epetra::Interface::Jacobian
{
	Teuchos::RCP<Epetra_Vector> oldX_;
	Teuchos::RCP<Epetra_Vector> B_;
	Teuchos::RCP<Epetra_Vector> xDot_;
	Teuchos::RCP<Ifpack_Preconditioner> precPtr_;

	double dt_;
	
public:
	SolverInterface(Epetra_Vector &newX, double dt)
		: dt_(dt)
		{
			oldX_    = Teuchos::rcp(new Epetra_Vector(newX));
			B_       = Teuchos::rcp(new Epetra_Vector(newX));
			xDot_    = Teuchos::rcp(new Epetra_Vector(newX));
		}
	
	~SolverInterface(){};             //destructor

	void setOldX(const Epetra_Vector &newX)	{*oldX_ = newX;}
	
	void setTimeStep(double dt) {dt_ = dt;}
	
	bool computeF(const Epetra_Vector &x, Epetra_Vector &f,
				  NOX::Epetra::Interface::Required::FillType F)
		{
			// Create non-owning rcp from Epetra_Vector f
			Teuchos::RCP<Epetra_Vector> f_rcp = Teuchos::rcp(&f, false);
			// Calculate rhs and store it in f
			THCM::Instance().evaluate(x, f_rcp, false);
			// Calculate mass matrix
			THCM::Instance().evaluateB();
			// Obtain mass matrix
			*B_ = THCM::Instance().DiagB();
			// Calculate d/dt x = (xnew - xold)/dt
			xDot_->Update(1.0/dt_, x, -1.0/dt_, *oldX_, 0.0);

			// Obtain number of local elements
			int NumMyElements     = xDot_->Map().NumMyElements();
			// Obtain mapping from local to global elements
			int *MyGlobalElements = xDot_->Map().MyGlobalElements();
			int index;
			double value;

			// Scale xdot with the mass matrix B
			for (int i = 0; i != NumMyElements; ++i)
			{
				index = MyGlobalElements[i];
				value = (*B_)[i] * (*xDot_)[i];
				xDot_->ReplaceGlobalValues(1, &value, MyGlobalElements + i);
			}
			// The final timestepping rhs is given by   B d/dt x + F(x) 
			f_rcp->Update(1.0, *xDot_, 1.0);
			return true;
		}	

	bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
		{
			THCM::Instance().evaluate(x, Teuchos::null, true);
			Epetra_CrsMatrix &JacCrs = dynamic_cast<Epetra_CrsMatrix &>(Jac);
			Teuchos::RCP<Epetra_CrsMatrix> Jacp = Teuchos::rcp(&JacCrs, false);
			Jacp = THCM::Instance().getJacobian();
			
			*B_ = THCM::Instance().DiagB();
			
			int NumMyElements     = Jacp->Map().NumMyElements();
			int *MyGlobalElements = Jacp->Map().MyGlobalElements();
			int index;
			double diag;

			for (int i = 0; i != NumMyElements; ++i)
			{
				diag = (*B_)[i] / dt_;
				Jacp->SumIntoGlobalValues(MyGlobalElements[i], 1,
										  &diag, MyGlobalElements + i);
			}
			Jacp->FillComplete();
			
			return true;
		}
    
	void exportMatrices()
		{
			*B_ = THCM::Instance().DiagB();
			Teuchos::RCP<Epetra_CrsMatrix> J = THCM::Instance().getJacobian();
			Teuchos::RCP<Epetra_CrsMatrix> B =
				Teuchos::rcp(new Epetra_CrsMatrix(Copy, J->RowMap(), 1));
			int NumMyElements     = B->Map().NumMyElements();
			int *MyGlobalElements = B->Map().MyGlobalElements();
			int index;
			double diag;
            
// TEUCHOS_ASSERT DIE MAPS VAN B_ EN B!!!!!
			for (int i = 0; i != NumMyElements; ++i)
			{
//				index = MyGlobalElements[i];
				diag  = (*B_)[i];
				B->InsertGlobalValues(MyGlobalElements[i], 1,
									  &diag, MyGlobalElements + i);
			}//setdiagvalues
			B->FillComplete();
			EpetraExt::RowMatrixToMatrixMarketFile("J.txt", *J);
			EpetraExt::RowMatrixToMatrixMarketFile("B.txt", *B);
		}
	bool computePrecMatrix(const Epetra_Vector &x, Epetra_RowMatrix &M)
		{
			std::cout << "*ERR* SolverInterface::computePrecMatrix()\n";
			std::cout << "*ERR* not using explicit preconditioning" << std::endl;
			exit(0);
			throw 1;
		}

	bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &O)
		{
			std::cout << "*ERR* SolverInterface::computePreconditioner()\n";
			std::cout << "*ERR* not using explicit preconditioning" << std::endl;
			exit(0);
			throw 1;
		}
};
#endif
