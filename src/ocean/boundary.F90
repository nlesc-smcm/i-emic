!*****************************************************************************
subroutine boundaries
  USE m_mat
  USE m_usr
  USE m_atm
  implicit none
  integer find_row2

  integer ii, i, j, k
  integer east, west, north, south, center
  integer neast,nwest,southw,southe, top, bottom
  integer eastb, westb, northb, southb
  integer neastb, nwestb, southwb, southeb
  integer eastt, westt, northt, southt
  integer neastt, nwestt, southwt, southet
  integer southee, easteast, northee, nnwest, nnorth, nneast
  integer nnorthee

  call TIMER_START('boundaries' // char(0))

  !  new stencil:
  !    +----------++-------++----------+
  !    | 12 15 18 || 3 6 9 || 21 24 27 |
  !    | 11 14 17 || 2 5 8 || 20 23 26 |
  !    | 10 13 16 || 1 4 7 || 19 22 25 |
  !    |  below   || center||  above   |
  !    +----------++-------++----------+

  ! Iterate over the flow domain
  do i = 1, n
     do j = 1, m
        do k = 1, l

           ! Give all the neighbours appropriate names.
           ! The landmask contains additional dummy cells on all borders.
           southw    = landm(i-1,j-1,k  )   !  1
           west      = landm(i-1,j  ,k  )   !  2
           nwest     = landm(i-1,j+1,k  )   !  3
           south     = landm(i  ,j-1,k  )   !  4
           center    = landm(i  ,j  ,k  )   !  5
           north     = landm(i  ,j+1,k  )   !  6
           southe    = landm(i+1,j-1,k  )   !  7
           east      = landm(i+1,j  ,k  )   !  8
           neast     = landm(i+1,j+1,k  )   !  9
           southwb   = landm(i-1,j-1,k-1)   ! 10
           westb     = landm(i-1,j  ,k-1)   ! 11
           nwestb    = landm(i-1,j+1,k-1)   ! 12
           southb    = landm(i  ,j-1,k-1)   ! 13
           bottom    = landm(i  ,j  ,k-1)   ! 14
           northb    = landm(i  ,j+1,k-1)   ! 15
           southeb   = landm(i+1,j-1,k-1)   ! 16
           eastb     = landm(i+1,j  ,k-1)   ! 17
           neastb    = landm(i+1,j+1,k-1)   ! 18
           southwt   = landm(i-1,j-1,k+1)   ! 19
           westt     = landm(i-1,j  ,k+1)   ! 20
           nwestt    = landm(i-1,j+1,k+1)   ! 21
           southt    = landm(i  ,j-1,k+1)   ! 22
           top       = landm(i  ,j  ,k+1)   ! 23
           northt    = landm(i  ,j+1,k+1)   ! 24
           southet   = landm(i+1,j-1,k+1)   ! 25
           eastt     = landm(i+1,j  ,k+1)   ! 26
           neastt    = landm(i+1,j+1,k+1)   ! 27

           if (i.lt.n) then

              ! Additional neighbours in the interior
              southee  = landm(i+2,j-1,k  ) !
              easteast = landm(i+2,j  ,k  ) !
              northee  = landm(i+2,j+1,k  ) !
              if (j.lt.m) then
                 nnorthee = landm(i+2,j+2,k)
              endif
           endif
           if (j.lt.m) then
              nnwest   = landm(i  ,j+2,k  ) !
              nnorth   = landm(i  ,j+2,k  ) !
              nneast   = landm(i  ,j+2,k  ) !
           endif

           !------- CENTER = OCEAN ---------------------------------------------------
           if (center == OCEAN) then
              ! check bottom cell first in order to prevent double checks
              ! on mirror/boundary points
              if (bottom == LAND) then  ! 14
                 if ((westb==LAND).and.(southwb==LAND).and.(southb==LAND)) then
                    An( 1,: ,UU,i,j,k) = An(1,: ,UU,i,j,k) + An(10,: ,UU,i,j,k) ! ACdN
                    An( 1,: ,VV,i,j,k) = An(1,: ,VV,i,j,k) + An(10,: ,VV,i,j,k) ! ACdN
                 endif
                 An(10,: ,UU,i,j,k) = 0.0
                 An(10,: ,VV,i,j,k) = 0.0
                 if ((westb==LAND).and.(neastb==LAND).and.(northb==LAND)) then
                    An( 2,: ,UU,i,j,k) = An(2,: ,UU,i,j,k) + An(11,: ,UU,i,j,k) ! ACdN
                    An( 2,: ,VV,i,j,k) = An(2,: ,VV,i,j,k) + An(11,: ,VV,i,j,k) ! ACdN
                 endif
                 An(11,: ,UU,i,j,k) = 0.0 !ACdN
                 An(11,: ,VV,i,j,k) = 0.0 !ACdN
                 if ((eastb==LAND).and.(southeb==LAND).and.(southb==LAND)) then
                    An( 4,: ,UU,i,j,k) = An(4,: ,UU,i,j,k) + An(13,: ,UU,i,j,k) ! ACdN
                    An( 4,: ,VV,i,j,k) = An(4,: ,VV,i,j,k) + An(13,: ,VV,i,j,k) ! ACdN
                 endif
                 An(13,: ,UU,i,j,k) = 0.0 !ACdN
                 An(13,: ,VV,i,j,k) = 0.0 !ACdN
                 if ((eastb==LAND).and.(neastb==LAND).and.(northb==LAND)) then
                    An( 5,: ,UU,i,j,k) = An(5,: ,UU,i,j,k) + An(14,: ,UU,i,j,k) ! ACdN
                    An( 5,: ,VV,i,j,k) = An(5,: ,VV,i,j,k) + An(14,: ,VV,i,j,k) ! ACdN
                 endif
                 An( 5,: ,TT,i,j,k) = An(5,: ,TT,i,j,k) + An(14,: ,TT,i,j,k) ! ACdN
                 An( 5,: ,SS,i,j,k) = An(5,: ,SS,i,j,k) + An(14,: ,SS,i,j,k) ! ACdN
                 An(14,: ,: ,i,j,k) = 0.0
              endif
              if (southwb == LAND) then ! 10
                 An(10,: ,: ,i,j,k) = 0.0
              endif
              if (westb == LAND) then   ! 11
                 An(11,: ,: ,i,j,k) = 0.0
              endif
              if (nwestb == LAND) then  ! 12
                 An(12,: ,: ,i,j,k) = 0.0
              endif
              if (southb == LAND) then  ! 13
                 An(13,: ,:,i,j,k)   = 0.0
              endif
              if (northb == LAND) then  ! 15
                 An(15,: ,:,i,j,k)  = 0.0
              endif
              if (southeb == LAND) then ! 16
                 An(16,: ,:,i,j,k)  = 0.0
              endif
              if (eastb == LAND) then   ! 17
                 An(17,: ,:,i,j,k)  = 0.0
              endif
              if (neastb == LAND) then  ! 18
                 An(18,: ,:,i,j,k)  = 0.0
              endif
              if (top == LAND) then ! 23
                 ! cannot occur in real flow domain, LAND above OCEAN is illegal
                 if (k.lt.(l)) then
                    write(f99,*) "NB error in boundary.f: LAND above OCEAN"
                    write(f99,*) i,j,k, landm(i  ,j  ,k),landm(i  ,j  ,k+1)
                 endif
                 if ((westt==LAND).and.(southwt==LAND).and.(southt==LAND)) then
                    An( 1,: ,UU,i,j,k) = An(1,: ,UU,i,j,k) + An(19,: ,UU,i,j,k) ! ACdN
                    An( 1,: ,VV,i,j,k) = An(1,: ,VV,i,j,k) + An(19,: ,VV,i,j,k) ! ACdN
                 endif
                 An(19,: ,UU,i,j,k) = 0.0
                 An(19,: ,VV,i,j,k) = 0.0
                 if ((westt==LAND).and.(nwestt==LAND).and.(northt==LAND)) then
                    An( 2,: ,UU,i,j,k) = An(2,: ,UU,i,j,k) + An(20,: ,UU,i,j,k) ! ACdN
                    An( 2,: ,VV,i,j,k) = An(2,: ,VV,i,j,k) + An(20,: ,VV,i,j,k) ! ACdN
                 endif
                 An(20,: ,UU,i,j,k) = 0.0 !ACdN
                 An(20,: ,VV,i,j,k) = 0.0 !ACdN
                 if ((eastt==LAND).and.(southet==LAND).and.(southt==LAND)) then
                    An( 4,: ,UU,i,j,k) = An(4,: ,UU,i,j,k) + An(22,: ,UU,i,j,k) ! ACdN
                    An( 4,: ,VV,i,j,k) = An(4,: ,VV,i,j,k) + An(22,: ,VV,i,j,k) ! ACdN
                 endif
                 An(22,: ,UU,i,j,k) = 0.0 !ACdN
                 An(22,: ,VV,i,j,k) = 0.0 !ACdN
                 if ((eastt==LAND).and.(neastt==LAND).and.(northt==LAND)) then
                    An( 5,: ,UU,i,j,k) = An(5,: ,UU,i,j,k) + An(23,: ,UU,i,j,k) ! ACdN
                    An( 5,: ,VV,i,j,k) = An(5,: ,VV,i,j,k) + An(23,: ,VV,i,j,k) ! ACdN
                 endif
                 An( 5,: ,TT,i,j,k) = An(5,: ,TT,i,j,k) + An(23,: ,TT,i,j,k) ! ACdN
                 An( 5,: ,SS,i,j,k) = An(5,: ,SS,i,j,k) + An(23,: ,SS,i,j,k) ! ACdN
                 An(23,: ,: ,i,j,k) = 0.0

                 Frc(find_row2(i,j,k,WW)) = 0.0

                 An( :,WW,: ,i,j,k) = 0.0
                 ! FIXME preconditioner breakdown if we remove the
                 ! connections, hence we try to maintain the
                 ! connection but make it inactive with 1e-10
                 An( 5, :,WW,i,j,k) = 1.0e-10 !MdT !TEM
                 An( 6, :,WW,i,j,k) = 1.0e-10 !MdT !TEM
                 An( 8, :,WW,i,j,k) = 1.0e-10 !MdT !TEM
                 An( 9, :,WW,i,j,k) = 1.0e-10 !MdT !TEM
                 An( 5,WW,WW,i,j,k) = 1.0

              endif
              !     if (southwt == LAND) then   ! 19
              if (southwt == LAND) then  ! 19
                 An(19,: ,:,i,j,k)  = 0.0
              endif
              !     if (westt == LAND) then     ! 20
              if (westt == LAND) then     ! 20
                 An(20,: ,: ,i,j,k) = 0.0
              endif
              if (nwestt == LAND) then   ! 21
                 An(21,:,:,i,j,k)    = 0.0
              endif
              if (southt == LAND) then   ! 22
                 An(22,: ,:,i,j,k)   = 0.0
              endif
              if (northt == LAND) then   ! 24
                 An(24,: ,:,i,j,k)  = 0.0
              endif
              if (southet == LAND) then ! 25
                 An(25,: ,:,i,j,k)  = 0.0
              endif
              if (eastt  == LAND) then     ! 26
                 An(26,: ,:,i,j,k)  = 0.0
              endif
              if (neastt == LAND) then   ! 27
                 An(27,: ,:,i,j,k)  = 0.0
              endif
              if (southw == LAND) then  ! 1
                 An( 1,: ,UU,i,j,k) = 0.0
                 An( 1,: ,VV,i,j,k) = 0.0
              endif
              if (west == LAND) then    ! 2
                 An( 5,: ,TT,i,j,k) = An(5,: ,TT,i,j,k) + An(2,: ,TT,i,j,k) ! ACdN
                 An( 5,: ,SS,i,j,k) = An(5,: ,SS,i,j,k) + An(2,: ,SS,i,j,k) ! ACdN
                 !              An( 5,TT,TT,i,j,k) = An(5,TT,TT,i,j,k) + An(2,TT,TT,i,j,k)
                 !              An( 5,SS,SS,i,j,k) = An(5,SS,SS,i,j,k) + An(2,SS,SS,i,j,k)
                 !              An( 5,TT,SS,i,j,k) = An(5,TT,SS,i,j,k) + An(2,TT,SS,i,j,k)
                 !              An( 5,SS,TT,i,j,k) = An(5,SS,TT,i,j,k) + An(2,SS,TT,i,j,k)
                 An( 2,: ,: ,i,j,k) = 0.0
                 An( 1,: ,UU,i,j,k) = 0.0
                 An( 1,: ,VV,i,j,k) = 0.0
              endif
              if (nwest == LAND) then   ! 3
                 An( 2,: ,UU,i,j,k) = 0.0
                 An( 2,: ,VV,i,j,k) = 0.0
                 An( 3,: ,UU,i,j,k) = 0.0
                 An( 3,: ,VV,i,j,k) = 0.0
              elseif (j.lt.m) then
                 if (nnwest == LAND) then
                    An( 3,: ,UU,i,j,k) = 0.0
                    An( 3,: ,VV,i,j,k) = 0.0
                 endif
              endif
              if (south == LAND) then   ! 4
                 An( 5,: ,SS,i,j,k) = An(5,: ,SS,i,j,k) + An(4,: ,SS,i,j,k) ! ACdN
                 An( 5,: ,TT,i,j,k) = An(5,: ,TT,i,j,k) + An(4,: ,TT,i,j,k) ! ACdN
                 !   An( 5,TT,TT,i,j,k) = An(5,TT,TT,i,j,k) + An(4,TT,TT,i,j,k)
                 !   An( 5,SS,SS,i,j,k) = An(5,SS,SS,i,j,k) + An(4,SS,SS,i,j,k)
                 !   An( 5,TT,SS,i,j,k) = An(5,TT,SS,i,j,k) + An(4,TT,SS,i,j,k)
                 !   An( 5,SS,TT,i,j,k) = An(5,SS,TT,i,j,k) + An(4,SS,TT,i,j,k)
                 An( 4,: ,: ,i,j,k) = 0.0
                 An( 1,: ,UU,i,j,k) = 0.0
                 An( 1,: ,VV,i,j,k) = 0.0
              endif
              if (north == LAND) then   ! 6
                 An( 2,: ,UU,i,j,k) = 0.0
                 An( 2,: ,VV,i,j,k) = 0.0
                 !
                 ! continuity
                 !
                 An( 2,PP,UU,i,j,k) = 0.0
                 An( 2,PP,VV,i,j,k) = 0.0
                 An( 5,PP,UU,i,j,k) = 0.0
                 An( 5,PP,VV,i,j,k) = 0.0
                 !
                 ! theta momentum
                 !
                 Frc(find_row2(i,j,k,VV)) = 0.0
                 An( :,VV,: ,i,j,k) = 0.0
                 An( 5,: ,VV,i,j,k) = 0.0 ! ACdN
                 An( 5,VV,VV,i,j,k) = 1.0
                 !
                 ! phi momentum
                 !
                 Frc(find_row2(i,j,k,UU)) = 0.0
                 An( :,UU, :,i,j,k) = 0.0
                 An( 5,: ,UU,i,j,k) = 0.0 ! ACdN
                 An( 5,UU,UU,i,j,k) = 1.0
                 !
                 ! tracers
                 !
                 An( 5,: ,SS,i,j,k) = An(5,: ,SS,i,j,k) + An(6,: ,SS,i,j,k) ! ACdN
                 An( 5,: ,TT,i,j,k) = An(5,: ,TT,i,j,k) + An(6,: ,TT,i,j,k) ! ACdN
                 !   An( 5,TT,TT,i,j,k) = An(5,TT,TT,i,j,k) + An(6,TT,TT,i,j,k)
                 !   An( 5,SS,SS,i,j,k) = An(5,SS,SS,i,j,k) + An(6,SS,SS,i,j,k)
                 !   An( 5,TT,SS,i,j,k) = An(5,TT,SS,i,j,k) + An(6,TT,SS,i,j,k)
                 !   An( 5,SS,TT,i,j,k) = An(5,SS,TT,i,j,k) + An(6,SS,TT,i,j,k)
                 An( 6,: ,: ,i,j,k) = 0.0
              elseif (j.lt.m) then
                 if (nnorth == LAND) then
                    An( 3,: ,UU,i,j,k) = 0.0
                    An( 3,: ,VV,i,j,k) = 0.0
                    An( 6,: ,UU,i,j,k) = 0.0
                    An( 6,: ,VV,i,j,k) = 0.0
                 endif
              endif
              if (southe == LAND) then  ! 7
                 An( 4, :,UU,i,j,k) = 0.0
                 An( 4, :,VV,i,j,k) = 0.0
                 An( 7, :,UU,i,j,k) = 0.0
                 An( 7, :,VV,i,j,k) = 0.0
              elseif (i.lt.n) then
                 if (southee == LAND) then
                    An( 7,: ,UU,i,j,k) = 0.0
                    An( 7,: ,VV,i,j,k) = 0.0
                 endif
              endif
              if (east == LAND) then    ! 8
                 An( 4,: ,UU,i,j,k) = 0.0
                 An( 4,: ,VV,i,j,k) = 0.0
                 !
                 ! continuity
                 !
                 An( 4,PP,UU,i,j,k) = 0.0
                 An( 4,PP,VV,i,j,k) = 0.0
                 An( 5,PP,UU,i,j,k) = 0.0
                 An( 5,PP,VV,i,j,k) = 0.0
                 !
                 ! phi momentum
                 !
                 Frc(find_row2(i,j,k,UU)) = 0.0
                 An( :,UU,: ,i,j,k) = 0.0
                 An( 5,: ,UU,i,j,k) = 0.0 ! ACdN
                 An( 5,UU,UU,i,j,k) = 1.0
                 !
                 ! theta momentum
                 !
                 Frc(find_row2(i,j,k,VV)) = 0.0
                 An( :,VV, :,i,j,k) = 0.0
                 An( 5,: ,VV,i,j,k) = 0.0 ! ACdN
                 An( 5,VV,VV,i,j,k) = 1.0
                 !
                 ! tracers
                 !
                 An( 5,: ,SS,i,j,k) = An(5,: ,SS,i,j,k) + An(8,: ,SS,i,j,k)
                 An( 5,: ,TT,i,j,k) = An(5,: ,TT,i,j,k) + An(8,: ,TT,i,j,k)
                 !   An( 5,TT,TT,i,j,k) = An(5,TT,TT,i,j,k) + An(8,TT,TT,i,j,k)
                 !   An( 5,SS,SS,i,j,k) = An(5,SS,SS,i,j,k) + An(8,SS,SS,i,j,k)
                 !   An( 5,TT,SS,i,j,k) = An(5,TT,SS,i,j,k) + An(8,TT,SS,i,j,k)
                 !   An( 5,SS,TT,i,j,k) = An(5,SS,TT,i,j,k) + An(8,SS,TT,i,j,k)
                 An( 8,: ,: ,i,j,k) = 0.0
                 An( 7, :,UU,i,j,k) = 0.0
                 An( 7, :,VV,i,j,k) = 0.0
              elseif (i.lt.n) then
                 if (easteast == LAND) then
                    An( 7,: ,UU,i,j,k) = 0.0
                    An( 7,: ,VV,i,j,k) = 0.0
                    An( 8,: ,UU,i,j,k) = 0.0
                    An( 8,: ,VV,i,j,k) = 0.0
                 endif
              endif
              if (neast == LAND) then   ! 9
                 !
                 ! phi momentum
                 !
                 Frc(find_row2(i,j,k,UU)) = 0.0
                 An( :,UU,: ,i,j,k) = 0.0
                 An( 5,: ,UU,i,j,k) = 0.0
                 An( 5,UU,UU,i,j,k) = 1.0
                 !
                 ! theta momentum
                 !
                 Frc(find_row2(i,j,k,VV)) = 0.0
                 An( :,VV,: ,i,j,k) = 0.0
                 An( 5,: ,VV,i,j,k) = 0.0
                 An( 5,VV,VV,i,j,k) = 1.0
                 An( 7, :,UU,i,j,k) = 0.0
                 An( 7, :,VV,i,j,k) = 0.0
              elseif ((i.lt.n).or.(j.lt.m)) then
                 if (i.lt.n) then
                    if (northee == LAND) then
                       An( 8,: ,UU,i,j,k) = 0.0
                       An( 8,: ,VV,i,j,k) = 0.0
                       An( 9,: ,UU,i,j,k) = 0.0
                       An( 9,: ,VV,i,j,k) = 0.0
                    elseif (j.lt.m) then
                       if (nnorthee == LAND) then
                          An( 9,: ,UU,i,j,k) = 0.0
                          An( 9,: ,VV,i,j,k) = 0.0
                       endif
                    endif
                 endif
                 if (j.lt.m) then
                    if (nneast == LAND) then
                       An( 6,: ,UU,i,j,k) = 0.0
                       An( 6,: ,VV,i,j,k) = 0.0
                       An( 9,: ,UU,i,j,k) = 0.0
                       An( 9,: ,VV,i,j,k) = 0.0
                    endif
                 endif
              endif
              
           else ! CENTER is not OCEAN so this should be on LAND
              An(:,:,:,i,j,k) = 0.0
              do ii = 1, nun
                 Frc(find_row2(i,j,k,ii)) = 0.0
                 An(5,ii,ii,i,j,k) = 1.0
              enddo
           endif
        enddo
     enddo
  enddo

  call TIMER_STOP('boundaries' // char(0))
end subroutine boundaries
