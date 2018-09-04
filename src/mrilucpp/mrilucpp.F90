!! interface to MRILU preconditioner, meant to be called from C++
module m_mrilucpp

#ifdef HAVE_IFPACK_MRILU
#ifndef WITH_UNION
#define csrmatrix anymatrix
#define prcmatrix anymatrix
#endif
  use m_build

  !! C-interoperability (Fortran 2003, I think...)
  use, intrinsic :: iso_c_binding !, only : c_double, c_int

  implicit none


  private
  public :: create,destroy,compute,apply,set_params

  !! double precision type
  integer, parameter :: dbl=8

  !! maximum number of instances
  integer, parameter :: max_inst = 25

  !! data type to hold:
  !! - the CRS matrix 
  !! - the MRILU preconditioner
  !! - all the MRILU parameters
  type :: MRILU_Prec

     type(csrmatrix),pointer :: A
     type(prcmatrix),pointer :: Prc

     integer :: blocksize

     LOGICAL ::  cutmck
     LOGICAL ::  scarow
     LOGICAL ::  xactelm
     LOGICAL ::  clsonce
     real(dbl)    ::  nlsfctr
     real(dbl)    ::  epsw
     real(dbl)    ::  elmfctr
     LOGICAL ::  gusmod
     real(dbl)    ::  gusfctr
     real(dbl)    ::  redfctr
     real(dbl)    ::  schtol
     real(dbl)    ::  denslim
     real(dbl)    ::  globfrac
     real(dbl)    ::  locfrac
     real(dbl)    ::  sparslim

     INTEGER ::  ilutype
     real(dbl)    ::  droptol
     real(dbl)    ::  compfct
     real(dbl)    ::  cpivtol
     real(dbl)    ::  lutol
     LOGICAL ::  singlu;
     INTEGER ::  outlev

     !! because of the way LOCA calls the construction routine
     !! we have to keep track of the allcoation status ourselves
     logical :: is_computed,is_created

  end type MRILU_Prec

  !! array of pointers to MRILU preconditioners instanciated from C++:
  type(mrilu_prec), dimension(max_inst) :: instance

  !! indicates which of the preconditioners are in use
  logical, dimension(max_inst) :: in_use=.false.

contains

  !! returns the first id which is not currently in use
  integer function first_unused()  

    implicit none

    integer :: i

    first_unused = -1

    do i=max_inst,1,-1
       if (.not. in_use(i)) then
          first_unused=i
       end if
    end do

  end function first_unused


  !! returns true if the index in question is valid (>0, <max_isnt) and in use
  logical function valid_id(id)

    implicit none

    integer :: id

    if ((id .le. 0) .or. (id .gt. max_inst))  then
       valid_id = .false.
    else
       valid_id = in_use(id)
    end if

  end function valid_id

  !! set the matrix to be ilu factored. It should be in C-style 
  !! (0-based) CRS format.                                      
  !! input:  id   - set this to 0 to create a new instance, or  
  !!                > 0 to overwrite an existing instance       
  !!         ndim - number of rows/cols of matrix               
  !!         nonz  - number of nonzero entries                  
  !!         beg  - CSR row pointer (size ndim+1)               
  !!         jco  - CSR column indices (size nnnz)              
  !!         jco  - CSR matrix entries (size nonz)              
  !!                                                            
  !! output: id, integer identifier of preconditioner instance, 
  !!         once assigned this should never be changed manually
  subroutine create(id,ndim,nonz,beg,jco,co) bind(C,name='mrilucpp_create')

    use m_wacsr

    implicit none

    integer(c_int), intent(inout) :: id 
    integer(c_int), intent(in)    :: ndim,nonz
    integer(c_int),dimension(ndim+1), intent(in) :: beg
    integer(c_int),dimension(nonz),intent(in)     :: jco
    real(c_double),dimension(nonz), intent(in)       :: co

    !_DEBUG2_('enter m_mriluprec::create, id=',id);

    if (id .eq. 0) then
       id = first_unused()
       if (id .lt. 0) then
          stop 'm_mriluprec::create: maximum number of MRILU instances exceeded!'
       end if
       in_use(id)=.true.
    else if (valid_id(id)) then
       call reset(id)
    else
       stop 'bad id passed into m_mriluprec::create'
    end if

    call wacsr(ndim,nonz,instance(id)%A)

    instance(id)%A%beg(1:ndim+1) = beg(1:ndim+1)+1
    instance(id)%A%jco(1:nonz) = jco(1:nonz)+1
    instance(id)%A%co(1:nonz) = co(1:nonz)

    instance(id)%is_created = .true.
    instance(id)%is_computed= .false.

    call default_params(id)

    !_DEBUG2_('leave m_mriluprec::create, id=',id);

  end subroutine create

  !! 'reset an instance for later use (destroy matrix and prec pointers)
  !! (private, not to be called from C++)
  subroutine reset(id) 

    use m_wfree
    implicit none

    integer(c_int), intent(in) :: id

    !_DEBUG2_('enter m_mriluprec::reset, id=',id);

    if (.not. valid_id(id)) then
       return
    end if

    if (instance(id)%is_computed) then
       call prcfree(instance(id)%Prc)
       ! note that the input matrix is reseted 
       ! by MRILU itself
    else if (instance(id)%is_created) then
       ! matrix was passed in but prec never computed.
       ! delete matrix
       call csrfree(instance(id)%A)
    end if

    instance(id)%is_created = .false.
    instance(id)%is_computed = .false.

    !_DEBUG2_('leave m_mriluprec::reset, id=',id);

  end subroutine reset

  !! destructor - reset precond and recycle its id
  !! (to be called from C++)
  subroutine destroy(id) bind(C,name='mrilucpp_destroy')

    implicit none

    integer(c_int) :: id

    if (valid_id(id)) then
       call reset(id)
       in_use(id) = .false.
    end if

    id=0 ! so the user can call create again

  end subroutine destroy

  !! convert integer to logical in a C++ consistent way
  logical function bool(int)

    implicit none

    integer, intent(in) :: int

    if (int==0) then
       bool=.false.
    else
       bool=.true.
    end if

  end function bool


  !! note: this horrible subroutine is not meant to be called by the 
  !! user but only by the C++ interface 'MRILU_Prec'
  subroutine set_params(id, blocksize_, cutmck_,scarow_,xactelm_,clsonce_,nlsfctr_,   &
       epsw_,elmfctr_,gusmod_,gusfctr_,redfctr_,schtol_,  &
       denslim_,globfrac_,locfrac_,sparslim_,ilutype_,    &
       droptol_,compfct_,cpivtol_,lutol_,singlu_,outlev_) &
       bind(C,name='mrilucpp_set_params')


    implicit none

    integer(c_int), intent(in) :: id
    real(c_double)    ::  nlsfctr_, epsw_, elmfctr_
    ! these are actually logicals, but as LOGICAL is not compatible with C++ bool,
    ! we pass them in as integers:
    INTEGER(c_int) ::  cutmck_, scarow_, xactelm_,clsonce_
    INTEGER(c_int) ::  gusmod_
    INTEGER(c_int) ::  singlu_

    real(c_double)    ::  gusfctr_, redfctr_, schtol_, denslim_, globfrac_
    real(c_double)    ::  locfrac_, sparslim_

    integer(c_int) ::  blocksize_, ilutype_
    real(c_double)    ::  droptol_, compfct_, cpivtol_, lutol_
    integer(c_int) ::  outlev_

    if (.not. valid_id(id)) then
       stop 'invalid id passed to m_mriluprec::set_params'
    end if

    call set_params_f(id, blocksize_, bool(cutmck_),bool(scarow_),bool(xactelm_),bool(clsonce_),&
         nlsfctr_, epsw_,elmfctr_,bool(gusmod_),gusfctr_,redfctr_,schtol_,  &
         denslim_,globfrac_,locfrac_,sparslim_,ilutype_,    &
         droptol_,compfct_,cpivtol_,lutol_,bool(singlu_),outlev_)

  end subroutine set_params

  !! this is the same as set_params, only with fortran-style logical args (do not use this
  !! from C++)
  subroutine set_params_f(id, blocksize_, cutmck_,scarow_,xactelm_,clsonce_,nlsfctr_,   &
       epsw_,elmfctr_,gusmod_,gusfctr_,redfctr_,schtol_,  &
       denslim_,globfrac_,locfrac_,sparslim_,ilutype_,    &
       droptol_,compfct_,cpivtol_,lutol_,singlu_,outlev_)


    implicit none

    integer, intent(in) :: id
    real(dbl)    ::  nlsfctr_, epsw_, elmfctr_
    LOGICAL ::  cutmck_, scarow_, xactelm_,clsonce_
    LOGICAL ::  gusmod_
    LOGICAL ::  singlu_

    real(dbl)    ::  gusfctr_, redfctr_, schtol_, denslim_, globfrac_
    real(dbl)    ::  locfrac_, sparslim_

    integer ::  blocksize_, ilutype_
    real(dbl)    ::  droptol_, compfct_, cpivtol_, lutol_
    integer ::  outlev_

    !_DEBUG2_('enter m_mriluprec::set_params_f, id=',id);

    if (.not. valid_id(id)) then
       stop 'invalid id passed to m_mriluprec::set_params'
    end if

    instance(id)%blocksize = blocksize_  ;
    instance(id)%outlev = outlev_  ;
    instance(id)%cutmck = cutmck_  ; instance(+id)%scarow =scarow_  ; 
    instance(id)%xactelm = xactelm_; instance(id)%clsonce=clsonce_  ;
    instance(id)%nlsfctr = nlsfctr_; instance(id)%epsw = epsw_      ; 
    instance(id)%elmfctr = elmfctr_; instance(id)%gusmod = gusmod_  ;
    instance(id)%gusfctr = gusfctr_; instance(id)%redfctr = redfctr_; 
    instance(id)%schtol = schtol_  ; instance(id)%denslim = denslim_;
    instance(id)%globfrac=globfrac_; instance(id)%locfrac = locfrac_; 
    instance(id)%sparslim=sparslim_; instance(id)%ilutype = ilutype_; 
    instance(id)%droptol = droptol_; instance(id)%compfct = compfct_; 
    instance(id)%cpivtol = cpivtol_; instance(id)%lutol = lutol_     ; 
    instance(id)%singlu = singlu_  ;           

    !_DEBUG2_('leave m_mriluprec::set_params_f, id=',id);

  end subroutine set_params_f

  !! set default parameters
  subroutine default_params(id)

    implicit none

    integer :: id

    INTEGER :: blocksize
    LOGICAL ::  cutmck, scarow, xactelm,clsonce
    real(dbl)    ::  nlsfctr, epsw, elmfctr
    LOGICAL ::  gusmod
    real(dbl)    ::  gusfctr, redfctr, schtol, denslim, globfrac
    real(dbl)    ::  locfrac, sparslim

    INTEGER ::  ilutype
    real(dbl)    ::  droptol, compfct, cpivtol, lutol
    LOGICAL ::  singlu
    INTEGER ::  outlev

    if (.not. valid_id(id)) then
       stop 'invalid id passed to m_mriluprec::default_params'
    end if

    blocksize = 1
    cutmck  = .false.
    scarow  = .false.
    xactelm = .false.
    clsonce = .false.
    nlsfctr = 1.0
    epsw    = 1.0e-2
    elmfctr = 0.1
    gusmod  = .false.
    gusfctr = 1.0
    redfctr = 0.8
    schtol  = 1.0e-6
    denslim = 1.0e-3
    globfrac= 1.0
    locfrac = 1.0
    sparslim= 0.675
    ilutype = 9
    droptol = 1.0e-4
    compfct = 1.0
    cpivtol = 0.1
    lutol   = 1.0e-14
    singlu  = .true.
    outlev  = 2

    call set_params_f(id, blocksize, cutmck,scarow,xactelm,clsonce,nlsfctr, &
         epsw,elmfctr,gusmod,gusfctr,redfctr,schtol,      &
         denslim,globfrac,locfrac,sparslim,ilutype,       &
         droptol,compfct,cpivtol,lutol,singlu,outlev)

  end subroutine default_params

  !! privat: pass params to MRILU module
  subroutine activate(id)

    use m_iniglb
    use m_iniprc   

    implicit none

    integer, intent(in) :: id

    ! the parameter setting routine is adopted from bgskit
    LOGICAL     :: a,b,c,d,t,u
    real(dbl)        :: e,f,g,h,i,j,k,l,m,n,p,q,r,s
    INTEGER     :: o

    if (.not. valid_id(id)) then
       stop 'invalid id passed to m_mriluprec::activate'
    end if

    a = instance(id)%cutmck;  b = instance(id)%scarow;   c = instance(id)%xactelm
    d = instance(id)%clsonce; e = instance(id)%nlsfctr;  f = instance(id)%epsw
    g = instance(id)%elmfctr; h = instance(id)%gusfctr;  i = instance(id)%redfctr;
    j = instance(id)%schtol;  k = instance(id)%denslim;  l = instance(id)%globfrac;
    m = instance(id)%locfrac; n = instance(id)%sparslim; o = instance(id)%ilutype;
    p = instance(id)%droptol; q = instance(id)%compfct;  r = instance(id)%cpivtol;
    s = instance(id)%lutol;   t = instance(id)%singlu;   u = .false.

#ifdef DEBUGGING
    write(*,*) "SETTING MRILU PARAMETERS: "
    write(*,*) "cutmck  = ",a
    write(*,*) "scarow  = ",b
    write(*,*) "xactelm = ",c
    write(*,*) "clsonce = ",d
    write(*,*) "nlsfctr = ",e
    write(*,*) "epsw    = ",f
    write(*,*) "elmfctr = ",g
    write(*,*) "gusfctr = ",h
    write(*,*) "redfctr = ",i
    write(*,*) "schtol  = ",j
    write(*,*) "denslim = ",k
    write(*,*) "globfrac = ",l
    write(*,*) "locfrac  = ",m
    write(*,*) "sparslim = ",n
    write(*,*) "ilutype  = ",o
    write(*,*) "droptol  = ",p
    write(*,*) "compfact = ",q
    write(*,*) "cpivtol  = ",r
    write(*,*) "lutol    = ",s
    write(*,*) "singlu   = ",t
    write(*,*) "TestPrec = ",u

#endif

    call iniglb(instance(id)%outlev)
    call iniprc (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)

  end subroutine activate

#ifdef DEBUGGING

  !! this subroutine is borrowed from Arie's 'sparsekit' module
  SUBROUTINE write_sparse_matrix(ixA,fname)
    !     write the matrix A to a file in asci format
    USE m_build
    IMPLICIT none
    !     INPUT/OUTPUT
    TYPE (csrmatrix), POINTER :: ixA
    character(len=*)          :: fname
    !     LOCAL
    INTEGER                   :: i,ier, n

    !_DEBUG_('enter m_mriluprec::write_sparse_matrix');
    !
    n = ixA%n
    !      OPEN(10,file = 'A.asc')
    !      WRITE(10,*) n
    !      DO i = 1, n+1
    !         WRITE(10,*) ixA%beg(i)
    !      END DO
    !      DO i = 1, ixA%beg(n+1)-1
    !         WRITE(10,*) ixA%jco(i)
    !      END DO
    !      DO i = 1, ixA%beg(n+1)-1
    !         WRITE(10,*) ixA%co(i)
    !      END DO
    !      CLOSE(10)
    OPEN(10,file=fname//'.info')
    OPEN(11,file=fname//'.beg')
    OPEN(12,file=fname//'.jco')
    OPEN(13,file=fname//'.co')
    WRITE(10,*) n,ixA%beg(n+1)-1
    DO i = 1,n+1
       write(11,*) ixA%beg(i)
    ENDDO
    DO i = 1,ixA%beg(n+1)-1
       write(12,*) ixA%jco(i)
    ENDDO
    DO i = 1,ixA%beg(n+1)-1
       write(13,*) ixA%co(i)
    ENDDO
    CLOSE(10)
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    !_DEBUG_('leave m_mriluprec::write_sparse_matrix');
  END SUBROUTINE write_sparse_matrix

#endif

  subroutine compute(id) bind(C,name='mrilucpp_compute')

    use m_build
    use m_wacsr
    use m_cmpprc
    use m_chkcnt

    implicit none

    integer(c_int), intent(in) :: id

    !_DEBUG2_('enter m_mriluprec::compute, id=',id);

    if (.not. valid_id(id)) then
       stop 'invalid id passed to m_mriluprec::compute'
    end if

    ! we do not have any reuse strategy here:
    if (.not. instance(id)%is_created) then
       write(*,*) "ERROR: no matrix available in mriluprec::compute!\n"
       write(*,*) "       'create' must be called first! \n"
       stop 'mriluprec error'
    end if

    if (instance(id)%is_computed) then
       write(*,*) 'ERROR: MRILU can be computed only once for a given matrix,'
       write(*,*) '       then the matrix has to be set again.               '
       stop 'mriluprec error'
    end if

#ifdef DEBUGGING
    call write_sparse_matrix(instance(id)%A,'MRILU_Matrix')
#endif

    ! set MRILU parameters...
    call activate(id)

    call cmpprc(instance(id)%blocksize,instance(id)%A,instance(id)%Prc)

    instance(id)%is_computed = .true.
    ! cmpprc resets the matrix:
    instance(id)%is_created = .false.

    !_DEBUG2_('leave m_mriluprec::compute, id=',id);

  end subroutine compute

  !! apply (inverse) preconditioner to a given vector vec.
  !! vec is overwritten with the solution.
  subroutine apply(id,ndim,rhs,sol) bind(C,name='mrilucpp_apply')

    use m_applprc

    implicit none

    integer(c_int), intent(in) :: id
    integer(c_int), intent(in) :: ndim
    real(c_double), dimension(ndim),intent(in) :: rhs
    real(c_double), dimension(ndim),intent(out) :: sol

    !_DEBUG2_('enter m_mriluprec::apply, id=',id);

    if (.not. valid_id(id)) then
       stop 'invalid id passed to m_mriluprec::apply'
    end if

    if (ndim .ne. instance(id)%Prc%n) then
       stop 'dimension mismatch in m_mriluprec::apply'
    end if

    if (.not. instance(id)%is_computed) then
       write(*,*) "WARNING: 'Apply' called before preconditioner was built!"
       write(*,*) "(",__FILE__,", line ",__LINE__,")"
       sol(1:ndim) = rhs(1:ndim)
    else
       sol(1:ndim) = 0.0
       CALL applprc(instance(id)%Prc,sol(1:ndim),rhs(1:ndim))
    end if

    !_DEBUG2_('leave m_mriluprec::apply, id=',id);

  end subroutine apply
#endif
end module m_mrilucpp
