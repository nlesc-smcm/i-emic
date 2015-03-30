
!! subroutines and data structures for monthly sea-surface temperature (SST),
!! salinity (SAL) and wind fields (TAU).
module m_monthly

use m_usr, usr_set_internal_forcing=>set_internal_forcing

implicit none


integer, parameter :: nt = 12 ! number of data sets (Jan-Dec)
real, dimension(:,:,:),allocatable :: mtaux,mtauy,mtatm,memip ! monthly mean fields
real, dimension(:,:,:,:),allocatable :: mtemp,msalt ! monthly mean fields
real, dimension(:,:),allocatable :: ataux,atauy,atatm,aemip ! annual mean fields
real, dimension(:,:,:),allocatable :: atemp,asalt  ! annual mean fields

logical :: initialized = .false.

contains

!! allocate arrays, read wind data. Levitus not implemented yet
subroutine init


implicit none

integer :: k

allocate(mtaux(n,m,nt),mtauy(n,m,nt),memip(n,m,nt),mtatm(n,m,nt));
allocate(mtemp(n,m,l,nt),msalt(n,m,l,nt));
allocate(ataux(n,m),atauy(n,m),aemip(n,m),atatm(n,m));
allocate(atemp(n,m,l),asalt(n,m,l));
write(f99,*) 'setup monthly forcing...'

ataux = taux
atauy = tauy
atatm = tatm
aemip = emip
atemp=internal_temp
asalt=internal_salt

! unless set_forcing is called, we use the annual mean for every month:

do k=1,nt
  mtaux(:,:,k) = ataux(:,:)
  mtauy(:,:,k) = atauy(:,:)
  mtatm(:,:,k) = atatm(:,:)
  memip(:,:,k) = aemip(:,:)
  mtemp(:,:,:,k)=atemp(:,:,:)
  msalt(:,:,:,k)=asalt(:,:,:)
  
end do

initialized = .true.

end subroutine init

!! set the levitus and wind field for a particular month.
!! The C++ code will loop over all months, create the    
!! fields in m_global, scatter them and pass them in     
!! through this subroutine.
subroutine set_forcing(ctatm,cemip,ctaux,ctauy,month)

use, intrinsic :: iso_c_binding
implicit none

real(c_double), dimension(m*n), intent(in) :: ctatm,cemip,ctaux,ctauy
integer(c_int), intent(in) :: month
integer :: i,j,pos

if (.not. initialized) then
  call init()
end if

if (month<=0 .or. month>nt) then
  stop 'ERROR in m_monthly::set_forcing: invalid month!'
end if

pos = 1
do j=1,m
  do i=1,n
    mtatm(i,j,month) = ctatm(pos)
    memip(i,j,month) = cemip(pos)
    mtaux(i,j,month) = ctaux(pos)
    mtauy(i,j,month) = ctauy(pos)
    pos = pos + 1
  end do
end do

end subroutine set_forcing

!! set the levitus and wind field for a particular month.
!! The C++ code will loop over all months, create the    
!! fields in m_global, scatter them and pass them in     
!! through this subroutine.
subroutine set_internal_forcing(ctemp,csalt,month)

use, intrinsic :: iso_c_binding
implicit none

real(c_double), dimension(m*n*l), intent(in) :: ctemp,csalt
integer(c_int), intent(in) :: month
integer :: i,j,k,pos

if (.not. initialized) then
  call init()
  end if
  
  if (month<=0 .or. month>nt) then
    stop 'ERROR in m_monthly::set_forcing: invalid month!'
    end if
    
    pos = 1
    do k=1,l
    do j=1,m
      do i=1,n
          mtemp(i,j,k,month) = ctemp(pos)
              msalt(i,j,k,month) = csalt(pos)
                  
                      pos = pos + 1
                        end do
                        end do
                        end do
                       
        end subroutine set_internal_forcing
        
!! update forcing (wind,surface salt/temperature) for   
!! time-dependent runs. Currently this only works for   
!! forcing from data (iza=ite=its=0). Data is available 
!! per month. The argument 't' is non-dimensional time. 
!! gamma is a continuation parameter, if the annual mean
!! field is called A and the monthly mean field M, then 
!! the field used is set to gamma*M + (1-gamma)*A.      
subroutine update_forcing(t, gammaW,gammaT,gammaS)

use, intrinsic :: iso_c_binding
use m_par
use m_usr
#if defined(INTEL)||defined(ORBIT)||defined(HOPF)
use ifport
#endif

implicit none

real(c_double) :: t
integer(c_int) :: i
integer(c_int) :: year
integer, dimension(4) :: month ! array indicating which months 
                               ! are to be used for interpolation
real, dimension(4) :: weight ! weights for interpolation
                              ! see split_time)

character(len=2) :: ibuf
character(len=1024) :: fname

real :: tatmmax,emipmax

real :: gammaW,gammaT,gammaS

! allocate data arrays if necessary (if this is done 
! here it means that set_forcing has not been called 
! and the monthly forcing is equal to the annual one)
if (.not. initialized) then
  call init()
end if

call split_time(t,year,month,weight)


! annual wind field component

taux(:,:) = (1.0-gammaW)*ataux
tauy(:,:) = (1.0-gammaW)*atauy

! annual temperature and salinity components
tatm(:,:) = (1.0-gammaT)*atatm
emip(:,:) = (1.0-gammaS)*aemip

! monthly components, interpolation coefficients are determined
!       by split_time routine
do i=1,4
  if (month(i)/=0) then
    taux(:,:) = taux(:,:)+gammaW*weight(i)*mtaux(:,:,month(i))
    tauy(:,:) = tauy(:,:)+gammaW*weight(i)*mtauy(:,:,month(i))
    tatm(:,:) = tatm(:,:)+gammaT*weight(i)*mtatm(:,:,month(i))
    emip(:,:) = emip(:,:)+gammaS*weight(i)*memip(:,:,month(i))
  end if
end do

end subroutine update_forcing

subroutine update_internal_forcing(t,gammaT,gammaS)

use, intrinsic :: iso_c_binding
use m_par
use m_usr
#if defined(INTEL)||defined(ORBIT)||defined(HOPF)
use ifport
#endif

implicit none

real(c_double) :: t
 real(c_double) :: gammaT,gammaS

integer :: i
integer :: year
integer, dimension(4) :: month ! array indicating which months 
                               ! are to be used for interpolation
  real, dimension(4) :: weight ! weights for interpolation
                                                             ! see split_time)
 character(len=2) :: ibuf
 character(len=1024) :: fname
 

                                                            
! annual temperature and salinity components
internal_temp(:,:,:) = (1.0-gammaT)*atemp
internal_salt(:,:,:) = (1.0-gammaS)*asalt
do i=1,4
  if (month(i)/=0) then
   
 internal_temp(:,:,:) =internal_temp(:,:,:) +gammaT*weight(i)*mtemp(:,:,:,month(i))
 internal_salt(:,:,:) =internal_salt(:,:,:) +gammaS*weight(i)*msalt(:,:,:,month(i))
   end if
     end do
    
  end subroutine update_internal_forcing

!! subroutine that takes non-dimensional time and
!! returns integers indicating the year and month
!! corresponding to that time value.            
!!                                              
!! to allow convenient interpolation, we now    
!! return 4 months to be used for the interpol.,
!! and the corresponding weights. The type of   
!! interpolation is set internally in this      
!! routine and has to be changed manually if    
!! desired.                                     
subroutine split_time(time,year,month,weight)

use m_par

implicit none

! here we choose the type of interpolation
integer, parameter :: itpl_type = 1     ! 0: constant
                                        ! 1: linear
                                        ! 2: 4-point (like in POP)
real, intent(in) :: time
integer,intent(out) :: year
integer, dimension(4) :: month
real, dimension(4) :: weight

real,parameter :: timesc=r0dim/udim

real :: time_in_secs
real :: secs_per_year, secs_per_month
integer :: this_month

secs_per_year  = 3600*24*365;
secs_per_month = secs_per_year/12.0

time_in_secs=time*timesc

! get time in years, months, days (int)
year = time_in_secs/secs_per_year
this_month = (time_in_secs - year*secs_per_year)/secs_per_month+1

month(:) = 0
weight(:) = 0.0

if (itpl_type==0) then ! piecewise constant interpolation
    month(1) = this_month
    weight(1) = 1.0
else if (itpl_type==1) then ! piecewise linear interpolation
  month(1) = this_month
  month(2) = mod(this_month+1,nt)
  if (month(2)==0) month(2)=nt
   weight(1) =  (time_in_secs-(year*secs_per_year+(this_month)*secs_per_month))/(-secs_per_month)
  weight(2) = 1.0-weight(1)
else if (itpl_type==2) then ! 4-point interpol, available in POP
  write(*,*) ' feature not implemented, ',__FILE__,', line ',__LINE__
  stop 'internal THCM error!'
else
  write(*,*) 'unknown itpl_type, ',__FILE__,', line ',__LINE__
  stop 'internal THCM error!'
end if
end subroutine split_time


end module m_monthly
