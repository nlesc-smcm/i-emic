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
        do k = 1, l+la

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
                    An(i,j,k, 1,: ,UU) = An(i,j,k,1,: ,UU) + An(i,j,k,10,: ,UU) ! ACdN
                    An(i,j,k, 1,: ,VV) = An(i,j,k,1,: ,VV) + An(i,j,k,10,: ,VV) ! ACdN
                 endif
                 An(i,j,k,10,: ,UU) = 0.0
                 An(i,j,k,10,: ,VV) = 0.0
                 if ((westb==LAND).and.(neastb==LAND).and.(northb==LAND)) then
                    An(i,j,k, 2,: ,UU) = An(i,j,k,2,: ,UU) + An(i,j,k,11,: ,UU) ! ACdN
                    An(i,j,k, 2,: ,VV) = An(i,j,k,2,: ,VV) + An(i,j,k,11,: ,VV) ! ACdN
                 endif
                 An(i,j,k,11,: ,UU) = 0.0 !ACdN
                 An(i,j,k,11,: ,VV) = 0.0 !ACdN
                 if ((eastb==LAND).and.(southeb==LAND).and.(southb==LAND)) then
                    An(i,j,k, 4,: ,UU) = An(i,j,k,4,: ,UU) + An(i,j,k,13,: ,UU) ! ACdN
                    An(i,j,k, 4,: ,VV) = An(i,j,k,4,: ,VV) + An(i,j,k,13,: ,VV) ! ACdN
                 endif
                 An(i,j,k,13,: ,UU) = 0.0 !ACdN
                 An(i,j,k,13,: ,VV) = 0.0 !ACdN
                 if ((eastb==LAND).and.(neastb==LAND).and.(northb==LAND)) then
                    An(i,j,k, 5,: ,UU) = An(i,j,k,5,: ,UU) + An(i,j,k,14,: ,UU) ! ACdN
                    An(i,j,k, 5,: ,VV) = An(i,j,k,5,: ,VV) + An(i,j,k,14,: ,VV) ! ACdN
                 endif
                 An(i,j,k, 5,: ,TT) = An(i,j,k,5,: ,TT) + An(i,j,k,14,: ,TT) ! ACdN
                 An(i,j,k, 5,: ,SS) = An(i,j,k,5,: ,SS) + An(i,j,k,14,: ,SS) ! ACdN
                 An(i,j,k,14,: ,: ) = 0.0
              endif
              if (southwb == LAND) then ! 10
                 An(i,j,k,10,: ,: ) = 0.0
              endif
              if (westb == LAND) then   ! 11
                 An(i,j,k,11,: ,: ) = 0.0
              endif
              if (nwestb == LAND) then  ! 12
                 An(i,j,k,12,: ,: ) = 0.0
              endif
              if (southb == LAND) then  ! 13
                 An(i,j,k,13,: ,:)   = 0.0
              endif
              if (northb == LAND) then  ! 15
                 An(i,j,k,15,: ,:)  = 0.0
              endif
              if (southeb == LAND) then ! 16
                 An(i,j,k,16,: ,:)  = 0.0
              endif
              if (eastb == LAND) then   ! 17
                 An(i,j,k,17,: ,:)  = 0.0
              endif
              if (neastb == LAND) then  ! 18
                 An(i,j,k,18,: ,:)  = 0.0
              endif
              if (top == LAND) then ! 23
                 ! cannot occur in real flow domain, LAND above OCEAN is illegal
                 if (k.lt.(l+la)) then
                    write(f99,*) "NB error in boundary.f: LAND above OCEAN"
                    write(f99,*) i,j,k, landm(i  ,j  ,k),landm(i  ,j  ,k+1)
                 endif
                 if ((westt==LAND).and.(southwt==LAND).and.(southt==LAND)) then
                    An(i,j,k, 1,: ,UU) = An(i,j,k,1,: ,UU) + An(i,j,k,19,: ,UU) ! ACdN
                    An(i,j,k, 1,: ,VV) = An(i,j,k,1,: ,VV) + An(i,j,k,19,: ,VV) ! ACdN
                 endif
                 An(i,j,k,19,: ,UU) = 0.0
                 An(i,j,k,19,: ,VV) = 0.0
                 if ((westt==LAND).and.(nwestt==LAND).and.(northt==LAND)) then
                    An(i,j,k, 2,: ,UU) = An(i,j,k,2,: ,UU) + An(i,j,k,20,: ,UU) ! ACdN
                    An(i,j,k, 2,: ,VV) = An(i,j,k,2,: ,VV) + An(i,j,k,20,: ,VV) ! ACdN
                 endif
                 An(i,j,k,20,: ,UU) = 0.0 !ACdN
                 An(i,j,k,20,: ,VV) = 0.0 !ACdN
                 if ((eastt==LAND).and.(southet==LAND).and.(southt==LAND)) then
                    An(i,j,k, 4,: ,UU) = An(i,j,k,4,: ,UU) + An(i,j,k,22,: ,UU) ! ACdN
                    An(i,j,k, 4,: ,VV) = An(i,j,k,4,: ,VV) + An(i,j,k,22,: ,VV) ! ACdN
                 endif
                 An(i,j,k,22,: ,UU) = 0.0 !ACdN
                 An(i,j,k,22,: ,VV) = 0.0 !ACdN
                 if ((eastt==LAND).and.(neastt==LAND).and.(northt==LAND)) then
                    An(i,j,k, 5,: ,UU) = An(i,j,k,5,: ,UU) + An(i,j,k,23,: ,UU) ! ACdN
                    An(i,j,k, 5,: ,VV) = An(i,j,k,5,: ,VV) + An(i,j,k,23,: ,VV) ! ACdN
                 endif
                 An(i,j,k, 5,: ,TT) = An(i,j,k,5,: ,TT) + An(i,j,k,23,: ,TT) ! ACdN
                 An(i,j,k, 5,: ,SS) = An(i,j,k,5,: ,SS) + An(i,j,k,23,: ,SS) ! ACdN
                 An(i,j,k,23,: ,: ) = 0.0

                 Frc(find_row2(i,j,k,WW)) = 0.0
                 An(i,j,k, :,WW,: ) = 0.0
                 ! FIXME preconditioner breakdown if we remove the
                 ! connections, hence we try to maintain the
                 ! connection but make it inactive with 1e-10
                 An(i,j,k, 5, :,WW) = 1.0e-10 !MdT !TEM
                 An(i,j,k, 6, :,WW) = 1.0e-10 !MdT !TEM
                 An(i,j,k, 8, :,WW) = 1.0e-10 !MdT !TEM
                 An(i,j,k, 9, :,WW) = 1.0e-10 !MdT !TEM
                 An(i,j,k, 5,WW,WW) = 1.0

              endif
              if (top == ATMOS) then ! deprecated in i-emic
                 An(i,j,k, 1,: ,UU) = An(i,j,k,1,: ,UU) + An(i,j,k,19,: ,UU) ! ACdN
                 An(i,j,k, 1,: ,VV) = An(i,j,k,1,: ,VV) + An(i,j,k,19,: ,VV) ! ACdN
                 An(i,j,k,19,: ,UU) = 0.0
                 An(i,j,k,19,: ,VV) = 0.0
                 An(i,j,k, 2,: ,UU) = An(i,j,k,2,: ,UU) + An(i,j,k,20,: ,UU) ! ACdN
                 An(i,j,k, 2,: ,VV) = An(i,j,k,2,: ,VV) + An(i,j,k,20,: ,VV) ! ACdN
                 An(i,j,k,20,: ,UU) = 0.0 !ACdN
                 An(i,j,k,20,: ,VV) = 0.0 !ACdN
                 An(i,j,k, 4,: ,UU) = An(i,j,k,4,: ,UU) + An(i,j,k,22,: ,UU) ! ACdN
                 An(i,j,k, 4,: ,VV) = An(i,j,k,4,: ,VV) + An(i,j,k,22,: ,VV) ! ACdN
                 An(i,j,k,22,: ,UU) = 0.0 !ACdN
                 An(i,j,k,22,: ,VV) = 0.0 !ACdN
                 An(i,j,k, 5,: ,UU) = An(i,j,k,5,: ,UU) + An(i,j,k,23,: ,UU) ! ACdN
                 An(i,j,k, 5,: ,VV) = An(i,j,k,5,: ,VV) + An(i,j,k,23,: ,VV) ! ACdN
                 An(i,j,k, 5,SS,SS) = An(i,j,k,5,SS,SS) + An(i,j,k,23,SS,SS)
                 An(i,j,k,23,SS,SS) = 0.0

                 Frc(find_row2(i,j,k,WW)) = 0.0
                 An(i,j,k, :,WW,: ) = 0.0
                 ! preconditioner breakdown if we do this:
                 ! An(i,j,k, 5, :,WW) = 0.0 ! 1.0e-10 !MdT
                 ! An(i,j,k, 6, :,WW) = 0.0 ! 1.0e-10 !MdT
                 ! An(i,j,k, 8, :,WW) = 0.0 ! 1.0e-10 !MdT
                 ! An(i,j,k, 9, :,WW) = 0.0 ! 1.0e-10 !MdT
                 An(i,j,k, 5,WW,WW) = 1.0

              endif
              !     if (southwt == LAND) then   ! 19
              if ((southwt == LAND).OR.(southwt ==ATMOS)) then  ! 19
                 An(i,j,k,19,: ,:)  = 0.0
              endif
              !     if (westt == LAND) then     ! 20
              if ((westt == LAND).OR.(westt == ATMOS)) then     ! 20
                 An(i,j,k,20,: ,: ) = 0.0
              endif
              if ((nwestt == LAND).OR.(nwestt == ATMOS)) then   ! 21
                 An(i,j,k,21,:,:)    = 0.0
              endif
              if ((southt == LAND).OR.(southt == ATMOS)) then   ! 22
                 An(i,j,k,22,: ,:)   = 0.0
              endif
              if ((northt == LAND).OR.(northt == ATMOS)) then   ! 24
                 An(i,j,k,24,: ,:)  = 0.0
              endif
              if ((southet == LAND).OR.(southet == ATMOS)) then ! 25
                 An(i,j,k,25,: ,:)  = 0.0
              endif
              if ((eastt == LAND).OR.(eastt == ATMOS)) then     ! 26
                 An(i,j,k,26,: ,:)  = 0.0
              endif
              if ((neastt == LAND).OR.(neastt == ATMOS)) then   ! 27
                 An(i,j,k,27,: ,:)  = 0.0
              endif
              if (southw == LAND) then  ! 1
                 An(i,j,k, 1,: ,UU) = 0.0
                 An(i,j,k, 1,: ,VV) = 0.0
              endif
              if (west == LAND) then    ! 2
                 An(i,j,k, 5,: ,TT) = An(i,j,k,5,: ,TT) + An(i,j,k,2,: ,TT) ! ACdN
                 An(i,j,k, 5,: ,SS) = An(i,j,k,5,: ,SS) + An(i,j,k,2,: ,SS) ! ACdN
                 !              An(i,j,k, 5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,2,TT,TT)
                 !              An(i,j,k, 5,SS,SS) = An(i,j,k,5,SS,SS) + An(i,j,k,2,SS,SS)
                 !              An(i,j,k, 5,TT,SS) = An(i,j,k,5,TT,SS) + An(i,j,k,2,TT,SS)
                 !              An(i,j,k, 5,SS,TT) = An(i,j,k,5,SS,TT) + An(i,j,k,2,SS,TT)
                 An(i,j,k, 2,: ,: ) = 0.0
                 An(i,j,k, 1,: ,UU) = 0.0
                 An(i,j,k, 1,: ,VV) = 0.0
              endif
              if (nwest == LAND) then   ! 3
                 An(i,j,k, 2,: ,UU) = 0.0
                 An(i,j,k, 2,: ,VV) = 0.0
                 An(i,j,k, 3,: ,UU) = 0.0
                 An(i,j,k, 3,: ,VV) = 0.0
              elseif (j.lt.m) then
                 if (nnwest == LAND) then
                    An(i,j,k, 3,: ,UU) = 0.0
                    An(i,j,k, 3,: ,VV) = 0.0
                 endif
              endif
              if (south == LAND) then   ! 4
                 An(i,j,k, 5,: ,SS) = An(i,j,k,5,: ,SS) + An(i,j,k,4,: ,SS) ! ACdN
                 An(i,j,k, 5,: ,TT) = An(i,j,k,5,: ,TT) + An(i,j,k,4,: ,TT) ! ACdN
                 !   An(i,j,k, 5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,4,TT,TT)
                 !   An(i,j,k, 5,SS,SS) = An(i,j,k,5,SS,SS) + An(i,j,k,4,SS,SS)
                 !   An(i,j,k, 5,TT,SS) = An(i,j,k,5,TT,SS) + An(i,j,k,4,TT,SS)
                 !   An(i,j,k, 5,SS,TT) = An(i,j,k,5,SS,TT) + An(i,j,k,4,SS,TT)
                 An(i,j,k, 4,: ,: ) = 0.0
                 An(i,j,k, 1,: ,UU) = 0.0
                 An(i,j,k, 1,: ,VV) = 0.0
              endif
              if (north == LAND) then   ! 6
                 An(i,j,k, 2,: ,UU) = 0.0
                 An(i,j,k, 2,: ,VV) = 0.0
                 !
                 ! continuity
                 !
                 An(i,j,k, 2,PP,UU) = 0.0
                 An(i,j,k, 2,PP,VV) = 0.0
                 An(i,j,k, 5,PP,UU) = 0.0
                 An(i,j,k, 5,PP,VV) = 0.0
                 !
                 ! theta momentum
                 !
                 Frc(find_row2(i,j,k,VV)) = 0.0
                 An(i,j,k, :,VV,: ) = 0.0
                 An(i,j,k, 5,: ,VV) = 0.0 ! ACdN
                 An(i,j,k, 5,VV,VV) = 1.0
                 !
                 ! phi momentum
                 !
                 Frc(find_row2(i,j,k,UU)) = 0.0
                 An(i,j,k, :,UU, :) = 0.0
                 An(i,j,k, 5,: ,UU) = 0.0 ! ACdN
                 An(i,j,k, 5,UU,UU) = 1.0
                 !
                 ! tracers
                 !
                 An(i,j,k, 5,: ,SS) = An(i,j,k,5,: ,SS) + An(i,j,k,6,: ,SS) ! ACdN
                 An(i,j,k, 5,: ,TT) = An(i,j,k,5,: ,TT) + An(i,j,k,6,: ,TT) ! ACdN
                 !   An(i,j,k, 5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,6,TT,TT)
                 !   An(i,j,k, 5,SS,SS) = An(i,j,k,5,SS,SS) + An(i,j,k,6,SS,SS)
                 !   An(i,j,k, 5,TT,SS) = An(i,j,k,5,TT,SS) + An(i,j,k,6,TT,SS)
                 !   An(i,j,k, 5,SS,TT) = An(i,j,k,5,SS,TT) + An(i,j,k,6,SS,TT)
                 An(i,j,k, 6,: ,: ) = 0.0
              elseif (j.lt.m) then
                 if (nnorth == LAND) then
                    An(i,j,k, 3,: ,UU) = 0.0
                    An(i,j,k, 3,: ,VV) = 0.0
                    An(i,j,k, 6,: ,UU) = 0.0
                    An(i,j,k, 6,: ,VV) = 0.0
                 endif
              endif
              if (southe == LAND) then  ! 7
                 An(i,j,k, 4, :,UU) = 0.0
                 An(i,j,k, 4, :,VV) = 0.0
                 An(i,j,k, 7, :,UU) = 0.0
                 An(i,j,k, 7, :,VV) = 0.0
              elseif (i.lt.n) then
                 if (southee == LAND) then
                    An(i,j,k, 7,: ,UU) = 0.0
                    An(i,j,k, 7,: ,VV) = 0.0
                 endif
              endif
              if (east == LAND) then    ! 8
                 An(i,j,k, 4,: ,UU) = 0.0
                 An(i,j,k, 4,: ,VV) = 0.0
                 !
                 ! continuity
                 !
                 An(i,j,k, 4,PP,UU) = 0.0
                 An(i,j,k, 4,PP,VV) = 0.0
                 An(i,j,k, 5,PP,UU) = 0.0
                 An(i,j,k, 5,PP,VV) = 0.0
                 !
                 ! phi momentum
                 !
                 Frc(find_row2(i,j,k,UU)) = 0.0
                 An(i,j,k, :,UU,: ) = 0.0
                 An(i,j,k, 5,: ,UU) = 0.0 ! ACdN
                 An(i,j,k, 5,UU,UU) = 1.0
                 !
                 ! theta momentum
                 !
                 Frc(find_row2(i,j,k,VV)) = 0.0
                 An(i,j,k, :,VV, :) = 0.0
                 An(i,j,k, 5,: ,VV) = 0.0 ! ACdN
                 An(i,j,k, 5,VV,VV) = 1.0
                 !
                 ! tracers
                 !
                 An(i,j,k, 5,: ,SS) = An(i,j,k,5,: ,SS) + An(i,j,k,8,: ,SS)
                 An(i,j,k, 5,: ,TT) = An(i,j,k,5,: ,TT) + An(i,j,k,8,: ,TT)
                 !   An(i,j,k, 5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,8,TT,TT)
                 !   An(i,j,k, 5,SS,SS) = An(i,j,k,5,SS,SS) + An(i,j,k,8,SS,SS)
                 !   An(i,j,k, 5,TT,SS) = An(i,j,k,5,TT,SS) + An(i,j,k,8,TT,SS)
                 !   An(i,j,k, 5,SS,TT) = An(i,j,k,5,SS,TT) + An(i,j,k,8,SS,TT)
                 An(i,j,k, 8,: ,: ) = 0.0
                 An(i,j,k, 7, :,UU) = 0.0
                 An(i,j,k, 7, :,VV) = 0.0
              elseif (i.lt.n) then
                 if (easteast == LAND) then
                    An(i,j,k, 7,: ,UU) = 0.0
                    An(i,j,k, 7,: ,VV) = 0.0
                    An(i,j,k, 8,: ,UU) = 0.0
                    An(i,j,k, 8,: ,VV) = 0.0
                 endif
              endif
              if (neast == LAND) then   ! 9
                 !
                 ! phi momentum
                 !
                 Frc(find_row2(i,j,k,UU)) = 0.0
                 An(i,j,k, :,UU,: ) = 0.0
                 An(i,j,k, 5,: ,UU) = 0.0
                 An(i,j,k, 5,UU,UU) = 1.0
                 !
                 ! theta momentum
                 !
                 Frc(find_row2(i,j,k,VV)) = 0.0
                 An(i,j,k, :,VV,: ) = 0.0
                 An(i,j,k, 5,: ,VV) = 0.0
                 An(i,j,k, 5,VV,VV) = 1.0
                 An(i,j,k, 7, :,UU) = 0.0
                 An(i,j,k, 7, :,VV) = 0.0
              elseif ((i.lt.n).or.(j.lt.m)) then
                 if (i.lt.n) then
                    if (northee == LAND) then
                       An(i,j,k, 8,: ,UU) = 0.0
                       An(i,j,k, 8,: ,VV) = 0.0
                       An(i,j,k, 9,: ,UU) = 0.0
                       An(i,j,k, 9,: ,VV) = 0.0
                    elseif (j.lt.m) then
                       if (nnorthee == LAND) then
                          An(i,j,k, 9,: ,UU) = 0.0
                          An(i,j,k, 9,: ,VV) = 0.0
                       endif
                    endif
                 endif
                 if (j.lt.m) then
                    if (nneast == LAND) then
                       An(i,j,k, 6,: ,UU) = 0.0
                       An(i,j,k, 6,: ,VV) = 0.0
                       An(i,j,k, 9,: ,UU) = 0.0
                       An(i,j,k, 9,: ,VV) = 0.0
                    endif
                 endif
              endif
              !------- CENTER = ATMOSPHERE ------------------------------------------
           else if (center == ATMOS) then
              write(*,*) 'center is atmos'
              if (west == LAND) then
                 An(i,j,k,5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,2,TT,TT)
                 An(i,j,k,2,TT,TT) = 0.0
              endif
              if (east == LAND) then
                 An(i,j,k,5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,8,TT,TT)
                 An(i,j,k,8,TT,TT) = 0.0
              endif
              if (north == LAND) then
                 An(i,j,k,5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,6,TT,TT)
                 An(i,j,k,6,TT,TT) = 0.0
              endif
              if (south == LAND) then
                 An(i,j,k,5,TT,TT) = An(i,j,k,5,TT,TT) + An(i,j,k,4,TT,TT)
                 An(i,j,k,4,TT,TT) = 0.0
              endif
              !------- CENTER = not OCEAN or ATMOSPHERE -----------------------------
              ! so this is probably on LAND
           else
              An(i,j,k,:,:,:) = 0.0
              do ii = 1, nun
                 Frc(find_row2(i,j,k,ii)) = 0.0
                 An(i,j,k,5,ii,ii) = 1.0
              enddo
           endif
        enddo
     enddo
  enddo
end subroutine boundaries
