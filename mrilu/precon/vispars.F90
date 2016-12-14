MODULE m_vispars

!-----------------------------------------------------------------------

!     Visualisation parameters.

!     These parameters control the visualisation of the matrices used
!     in the construction of the Multi Level Preconditioner

!     The parameters get a default value
!     during declaration. The values of these
!     parameters can be changed by the SUBROUTINE 'inivis'.


!     Description of the parameters in /vispars/:
!     ===========================================
!     NROwrtn  Number ReOrdered matrices written to files.
!     NSCwrtn  Number Schur-complement matrices written to files.
!     VisNRO   Visualise the first 'visNRO' ReOrdered matrices.
!     VisNSC   Visualise the first 'visNSC' Schur-complement matrices.
!     VisAsc   Visualise the original matrix A (possibly scaled).
!     VisLSC   Visualise the last Schur-complement.
!     VisILDU  Visualize the Incomplete LDU factorization of the last
!              Schur-complement.

INTEGER ::			&
	nrowrtn	=0		,&
	nscwrtn	=0		,&
	visnro	=0		,& 
	visnsc	=0

LOGICAL :: 			&
	visasc	=.false.	,&
	visildu	=.false.	,&
	vislsc	=.false.


!-----------------------------------------------------------------------

END MODULE