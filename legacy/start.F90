	MODULE m_start
!          INCLUDE 'usr.com'
          !     COMMON
          logical :: fix
          integer :: cstep,NM,Nopt,PM
          real    :: Kn,Knc,Tn,Tnc
          real    :: ds,maxds,timea,timeb
          REAL, ALLOCATABLE, DIMENSION(:) :: rb2old, r2o, dfdmu, rh0,rh1,rh2
!          common /SP1/  rb2old(ndim),r2o(ndim),dfdmu(ndim),rh0(ndim),rh1(ndim),rh2(ndim)
!          common /SP2/  fix,cstep,timea,timeb
!          common /SP3/  Kn,Knc,Tn,Tnc,ds,maxds
        END MODULE m_start
