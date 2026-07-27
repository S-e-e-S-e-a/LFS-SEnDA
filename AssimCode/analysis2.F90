subroutine analysis2(A4,psi4,D4, R4, S4, ndim, nz, nrens, nrobs, alpha, assim_mask, mlon, mlat, anal_solver, vert_local, dep_obs_ind, type_obs)
  !lzt 20120312 add alpha
   use m_multa
   use m_chol
   use m_vert_loc 
   implicit none
   logical, intent(in) :: vert_local
   integer, intent(in) :: ndim, nz             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   integer, intent(in) :: dep_obs_ind(nrobs)
   character(len=3), intent(in) :: type_obs(nrobs)
   real, intent(inout) :: A4(ndim,nrens)   ! ensemble matrix
   real, intent(inout) :: psi4(ndim)       ! state vector (forecast -> vector)
   real, intent(in)    :: D4(nrobs)  ! matrix  holding observation innovations
   real, intent(in)    :: S4(nrobs,nrens)  ! matrix holding HA' 
   real, intent(inout) :: R4(nrobs,nrobs)  ! Error covariance matrix for observations
!   logical, intent(in) :: verbose
   real, intent(in):: alpha, assim_mask, mlon, mlat

!   real, allocatable, dimension(:,:) :: X1,X2,Reps,I_N,U,V
!   real*8, allocatable, dimension(:) :: X4,X5
   real*8, allocatable, dimension(:,:) :: X1,X2,I_N,U,V,Rinv
   real*8, allocatable, dimension(:,:) :: X3,X4,X5
   real*8, allocatable, dimension(:,:) :: A_IN,B,Lz,Lc
   real*8, allocatable, dimension(:)   :: sig,work

   real*8 :: A(ndim,nrens)
   real*8 :: D(nrobs,1),psi(ndim,1),d_psi(ndim,1) ,beta(ndim,1)
   real*8 :: S(nrobs,nrens)
   real*8 :: R(nrobs,nrobs)
   real*8 :: T_surf, weight
   real*8, parameter :: R_max = 100.0
   real*8, parameter :: k_steep = 5.0
   real*8, parameter :: T_th_NH = -0.5
   real*8, parameter :: T_th_SH = 1.5
   real*8, parameter :: pi = 3.14159265
   character(len=3) :: anal_solver
!   real*4 :: levels(nz)

   real*8 sigsum,sigsum1,oneobs(1,1)
   integer ierr,nrsigma,i,j,lwork,m,ii,jj,k
   integer iblkmax, iens,ipiv(nrobs),info1,info2
   integer, parameter :: target= 2*30+2
   character(len=2) tag2
   real :: h_weight, v_weight, damp
           A=DBLE(A4)
           psi(:,1)=DBLE(psi4)
           D(1:nrobs,1)=DBLE(D4(1:nrobs))
           S=DBLE(S4)
           R=DBLE(R4)
           R=R/alpha

    T_surf = psi(2*nz+2,1)
    weight = 1.0
    if (assim_mask == 0.0) then
       weight = 1.0 / R_max
    else if (mlat > 40.0) then
       weight = 1.0 / (1.0 + (R_max - 1.0) / (1.0 + exp(k_steep * (T_surf - T_th_NH))))
    else if (mlat < -40.0) then
       weight = 1.0 / (1.0 + (R_max - 1.0) / (1.0 + exp(k_steep * (T_surf - T_th_SH))))
    endif
    R = R / weight
    if (mlat >= 60.0) then
       weight = weight * 0.5 * (1.0 + cos(pi * min((mlat - 60.0)/15.0, 1.0)))
    else if (mlat <= -60.0) then
       weight = weight * 0.5 * (1.0 + cos(pi * min((-60.0 - mlat)/15.0, 1.0)))
    end if
    D = D * weight


!   if ((assim_mask == 0.0) .or. &
!      (psi(2*nz+2,1) < -1.5 .and. mlat  > 40.) .or. &
!         (psi(2*nz+2,1) <  1.2 .and. mlat < -40.)) then
!      R = 1600.0 * R
!    else if ((psi(2*nz+2,1) < -1.0 .and. mlat  > 40.) .or. &
!        (psi(2*nz+2,1) <  1.5 .and. mlat  < -40.)) then
!      R = 900.0 * R
!    else if ((psi(2*nz+2,1) < -0.5 .and. mlat  > 40.) .or. &
!        (psi(2*nz+2,1) <  1.8 .and. mlat  < -40.)) then
!      R = 400.0 * R
!    else if ((psi(2*nz+2,1) < 0 .and. mlat  > 40.) .or. &
!        (psi(2*nz+2,1) <  2. .and. mlat  < -40.)) then
!      R = 100.0 * R
!    endif

!   if (mlat>=72.) then
!        D = D*max((1.0-(mlat-72.0)/15.0), 0.0) ! D=0 in 87 N
!   else if (mlat<=-62.) then
!        D = D*max(((77.0+mlat)/15.0), 0.0)  ! D=0 in 77 S
!   end if

   !call get_nz(type_obs, lon_obs, lat_obs, nrobs, nz_obs)
   !do j=1,nrobs
   ! R(j,j) = R(j,j)*nz_obs(j)
    !S(j,:) = S(j,:)/dble(sqrt(nz_obs(j)))
    !D(j,1) = D(j,1)/dble(sqrt(nz_obs(j)))
   !end do   

   if (nrobs > 1) then
      R=dble(nrens)*R+matmul(S,transpose(S))
      if (sum(R) /= sum(R) .or. sum(R) > 1.0d30 .or. sum(R) < -1.0d30) then
         d_psi = 0.0d0
         goto 999
      end if
      if ( anal_solver == 'SVD' ) then
      allocate (U(nrobs,nrobs)  )
      allocate (V(nrobs,nrobs)  )
      allocate (sig(nrobs)  )
      lwork=2*max(3*nrobs+nrobs,5*nrobs)
      allocate(work(lwork))
      sig=0.0
      call dgesvd('A', 'A', nrobs, nrobs, R, nrobs, sig, U, nrobs, V, nrobs, work, lwork, ierr)
      deallocate(work)

      if (ierr /= 0) then
         print *,'ierr from call dgesvd= ',ierr
         stop
      endif

      sigsum=sum( sig(1:nrobs) )
      sigsum1=0.0
      nrsigma=0
      do i=1,min(nrobs,nrens)                 ! singular values are in descending order
         if (sigsum1/sigsum < 0.998) then 
            nrsigma=nrsigma+1
            sigsum1=sigsum1+sig(i)
            sig(i) = 1.0/sig(i)
         else
            sig(i:nrobs)=0.0
            exit
         endif
      enddo

      allocate (X1(nrsigma,nrobs))
      do i=1,nrsigma
      do j=1,nrobs
         X1(i,j)=sig(i)*U(j,i)
      enddo
      enddo
      deallocate(sig)
      allocate (X2(nrsigma,1))
      X2= matmul(X1,D)
      deallocate(X1) 
      allocate (X3(nrobs,1))
      X3=matmul(transpose(V(1:nrsigma,:)),X2)
      deallocate(U)
      deallocate(V)
      deallocate(X2)
      allocate (I_N(nrens,nrens))
      I_N=-1.0/dble(nrens)
      do i=1,nrens
         I_N(i,i)=I_N(i,i)+1.0
      enddo
      allocate(X4(nrens,1))
      X4=matmul(transpose(S),X3)
      deallocate(X3)
      allocate(X5(nrens,1))
      X5=matmul(I_N,X4)
      deallocate(X4)
      deallocate(I_N)
      d_psi=matmul(A,X5)
      deallocate(X5)
    
    else if ( anal_solver == 'CHO') then
      allocate (Rinv(nrobs,nrobs))
      call chol_inv(R, Rinv, nrobs)
      allocate(X1(nrobs,1))
      X1=matmul(Rinv,D)
      allocate(X2(nrens,1))
      X2=matmul(transpose(S),X1)
      deallocate(X1)
      d_psi=matmul(A,X2)
      deallocate(X2)
    else 
      d_psi=0
    endif
   else   ! ONLY 1 OBS
      oneobs=matmul(S,transpose(S))+R*dble(nrens)
      print *,'oneobs: ',oneobs(1,1)
      allocate (X3(nrobs,1))
      X3(1,1)=D(1,1)/oneobs(1,1)
      allocate (I_N(nrens,nrens))
      I_N=-1.0/dble(nrens)
      do i=1,nrens
         I_N(i,i)=I_N(i,i)+1.0
      enddo
      allocate(X4(nrens,1))
      X4=matmul(transpose(S),X3)
      deallocate(X3)
      allocate(X5(nrens,1))
      X5=matmul(I_N,X4)
      deallocate(X4)
      deallocate(I_N)
      d_psi=matmul(A,X5)
      deallocate(X5)
   endif
999 continue
    ! z0 uu vv [-2, 2]
    where (d_psi(1:2*nz+1,1) >  2.0) d_psi(1:2*nz+1,1) =  2.0
    where (d_psi(1:2*nz+1,1) < -2.0) d_psi(1:2*nz+1,1) = -2.0

    ! tt [-5, 5]
    where (d_psi(2*nz+2:3*nz+1,1) >  5.0) d_psi(2*nz+2:3*nz+1,1) =  5.0
    where (d_psi(2*nz+2:3*nz+1,1) < -5.0) d_psi(2*nz+2:3*nz+1,1) = -5.0

    ! ss [-6, 6]
    where (d_psi(3*nz+2:4*nz+1,1) >  6.0) d_psi(3*nz+2:4*nz+1,1) =  6.0
    where (d_psi(3*nz+2:4*nz+1,1) < -6.0) d_psi(3*nz+2:4*nz+1,1) = -6.0
        
!    if ((assim_mask == 0.0) .or. (T_surf < -1.5 .and. mlat > 40.) .or. &
!         (T_surf <  1.0 .and. mlat < -40.)) then
!      d_psi=0
!    end if
    
!    if ((mlon-257.75)**2.+(mlat-5.375)**2. < 50) then
!       h_weight = 1.0 - ((mlon-257.75)**2.+(mlat-5.375)**2.) / 50
!       do k = 1, 55
!          v_weight = real(min(k,5) - 1) / 4.0
!          damp = 1.0 - h_weight * v_weight
!          ! tt increment
!          d_psi(2*nz+1+k,1) = d_psi(2*nz+1+k,1) * damp
!          ! ss increment
!          d_psi(3*nz+1+k,1) = d_psi(3*nz+1+k,1) * damp
!       end do
!    end if
!    if ((mlon-240.5)**2.+(mlat+1.625)**2. < 18) then
!       h_weight = 1.0 - ((mlon-240.5)**2.+(mlat+1.625)**2.) / 18
!       do k = 1, 55
!          v_weight = real(min(k,5) - 1) / 4.0
!          damp = 1.0 - h_weight * v_weight
!          ! tt increment
!          d_psi(2*nz+1+k,1) = d_psi(2*nz+1+k,1) * damp
!          ! ss increment
!          d_psi(3*nz+1+k,1) = d_psi(3*nz+1+k,1) * damp
!       end do
!    end if
!    if ((mlon-309.0)**2.+(mlat+40.375)**2. < 50) then
!       h_weight = 1.0 - ((mlon-309.0)**2.+(mlat+40.375)**2.) / 50
!       do k = 1, 55
!          v_weight = real(min(k,5) - 1) / 4.0
!          damp = 1.0 - h_weight * v_weight
!          ! tt increment
!          d_psi(2*nz+1+k,1) = d_psi(2*nz+1+k,1) * damp
!          ! ss increment
!          d_psi(3*nz+1+k,1) = d_psi(3*nz+1+k,1) * damp
!       end do
!    end if

    psi=d_psi+psi

!lzt
    psi4=REAL(psi(:,1))
    A4=REAL(A)
    R4=REAL(R) 
     
end subroutine analysis2

subroutine analysis2_letkf(A4,psi4,D4, R4, S4, ndim, nz, nrens, nrobs, alpha, assim_mask, mlon, mlat, anal_solver, vert_local, dep_obs_ind, type_obs)
  !lzt 20120312 add alpha
   use m_multa
   use m_chol
   use m_vert_loc
   implicit none
   logical, intent(in) :: vert_local
   integer, intent(in) :: ndim, nz             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   integer, intent(in) :: dep_obs_ind(nrobs)
   character(len=3), intent(in) :: type_obs(nrobs)
   real, intent(inout) :: A4(ndim,nrens)   ! ensemble matrix
   real, intent(inout) :: psi4(ndim)       ! state vector (forecast -> vector)
   real, intent(in)    :: D4(nrobs)  ! matrix  holding observation innovations
   real, intent(in)    :: S4(nrobs,nrens)  ! matrix holding HA'
   real, intent(inout) :: R4(nrobs,nrobs)  ! Error covariance matrix for observations
!   logical, intent(in) :: verbose
   real, intent(in):: alpha, assim_mask, mlon, mlat
   real*8, allocatable, dimension(:,:) :: X1,X2,I_N,C,Cinv,Rinv,U,V
   real*8, allocatable, dimension(:,:) :: X3,X4,X5
   real*8, allocatable, dimension(:,:) :: A_IN,B,Lz,Lc
   real*8, allocatable, dimension(:)   :: sig,work
   real*8 :: A(ndim,nrens)
   real*8 :: D(nrobs,1),psi(ndim,1),d_psi(ndim,1) ,beta(ndim,1)
   real*8 :: S(nrobs,nrens)
   real*8 :: R(nrobs,nrobs)
   real*8 :: T_surf, weight
   real*8, parameter :: R_max = 100.0
   real*8, parameter :: k_steep = 5.0 
   real*8, parameter :: T_th_NH = -0.5
   real*8, parameter :: T_th_SH = 1.5  
   real*8, parameter :: pi = 3.14159265
   character(len=3) :: anal_solver
!   real*4 :: levels(nz)
!lzt--------------

   real*8 traceC,eps,lambda_mean,sigsum,sigsum1,oneobs(1,1)
   integer ierr,nrsigma,i,j,lwork,m,ii,jj,k
   integer iblkmax, iens,ipiv(nrobs),info1,info2
   integer, parameter :: target= 2*30+2  ! TEM(01) after h0,u,v
   character(len=2) tag2
!lzt
           A=DBLE(A4)
           psi(:,1)=DBLE(psi4)
           D(1:nrobs,1)=DBLE(D4(1:nrobs))
           S=DBLE(S4)
           R=DBLE(R4)
    
    T_surf = psi(2*nz+2,1)
    weight = 1.0
    if (assim_mask == 0.0) then
       weight = 1.0 / R_max
    else if (mlat > 40.0) then
       weight = 1.0 / (1.0 + (R_max - 1.0) / (1.0 + exp(k_steep * (T_surf - T_th_NH))))
    else if (mlat < -40.0) then
       weight = 1.0 / (1.0 + (R_max - 1.0) / (1.0 + exp(k_steep * (T_surf - T_th_SH))))
    endif
    R = R / weight
    if (mlat >= 60.0) then
       weight = weight * 0.5 * (1.0 + cos(pi * min((mlat - 60.0)/15.0, 1.0)))
    else if (mlat <= -60.0) then
       weight = weight * 0.5 * (1.0 + cos(pi * min((-60.0 - mlat)/15.0, 1.0)))
    end if
    D = D * weight

!   if ((assim_mask == 0.0) .or. &
!      (psi(2*nz+2,1) < -1.5 .and. mlat  > 40.) .or. &
!         (psi(2*nz+2,1) <  1.2 .and. mlat < -40.)) then
!      R = 1600.0 * R
!    else if ((psi(2*nz+2,1) < -1.0 .and. mlat  > 40.) .or. &
!        (psi(2*nz+2,1) <  1.5 .and. mlat  < -40.)) then
!      R = 900.0 * R
!    else if ((psi(2*nz+2,1) < -0.5 .and. mlat  > 40.) .or. &
!        (psi(2*nz+2,1) <  1.8 .and. mlat  < -40.)) then
!      R = 400.0 * R
!    else if ((psi(2*nz+2,1) < 0 .and. mlat  > 40.) .or. &
!        (psi(2*nz+2,1) <  2. .and. mlat  < -40.)) then
!      R = 100.0 * R
!    endif

!   if (mlat>=72.) then
!        D = D*max((1.0-(mlat-72.0)/15.0), 0.0) ! D=0 in 87 N
!   else if (mlat<=-62.) then
!        D = D*max(((77.0+mlat)/15.0), 0.0)  ! D=0 in 77 S
!   end if

   do i = 1, nrens
      A(:, i) = A(:, i) - sum(A, dim=2) / dble(nrens)
      S(:, i) = S(:, i) - sum(S, dim=2) / dble(nrens)
   end do
   
   if (nrobs > 1) then
      allocate (I_N(nrens,nrens))
      allocate (C(nrens,nrens))
      allocate (Rinv(nrobs,nrobs))
      I_N = 0.0d0
      do i = 1, nrens
        I_N(i,i) = 1.0d0
      end do
      Rinv = 0.0d0 
      do i = 1, nrobs
         Rinv(i,i) = 1.0d0 / max(R(i,i), 1.0d-8)
      end do
      C=alpha*(nrens-1)*I_N+matmul(matmul(transpose(S),Rinv),S)
      deallocate(I_N)
      if (sum(C) /= sum(C) .or. sum(C) > 1.0d30 .or. sum(C) < -1.0d30) then
         d_psi = 0.0d0
         deallocate(C)
         deallocate(Rinv)
         goto 999
      end if
      if ( anal_solver == 'SVD' ) then
      allocate (U(nrens,nrens)  )
      allocate (V(nrens,nrens)  )
      allocate (sig(nrens)  )
      lwork=2*max(3*nrens+nrens,5*nrens)
      allocate(work(lwork))
      sig=0.0
      call dgesvd('A', 'A', nrens, nrens, C, nrens, sig, U, nrens, V, nrens, work, lwork, ierr)
      deallocate(work)

      if (ierr /= 0) then
         print *,'ierr from call dgesvd= ',ierr
         stop
      endif

      sigsum=sum( sig(1:nrens) )
      sigsum1=0.0
      nrsigma=0
      do i=1,nrens             ! singular values are in descending order
         if (sigsum1/sigsum < 0.999) then
            nrsigma=nrsigma+1
            sigsum1=sigsum1+sig(i)
            sig(i) = 1.0/sig(i)
         else
            sig(i:nrens)=0.0
            exit
         endif
      enddo
      deallocate(C)
      allocate (X1(nrsigma,nrens))
      do i=1,nrsigma
      do j=1,nrens
         X1(i,j)=sig(i)*U(j,i)
      enddo
      enddo
      deallocate(sig)
      allocate (X2(nrobs,1))
      X2= matmul(Rinv,D)
      allocate (X3(nrens,1))
      X3= matmul(transpose(S),X2)
      deallocate(X2)
      allocate (X4(nrsigma,1))
      X4=matmul(X1,X3)
      deallocate(X1,X3)
      allocate (X5(nrens,1))
      X5=matmul(transpose(V(1:nrsigma,:)),X4)
      deallocate(U,V,X4)
      d_psi=matmul(A,X5)
      deallocate(X5)

    else if ( anal_solver == 'CHO') then
      allocate(Cinv(nrens,nrens))
      call chol_inv(C, Cinv, nrens)
      allocate (X2(nrobs,1))
      X2=matmul(Rinv,D)
      deallocate(C)
      deallocate(Rinv)
      allocate (X3(nrens,1))
      X3=matmul(transpose(S),X2)
      allocate (X4(nrens,1))
      X4=matmul(Cinv,X3)
      deallocate(Cinv)
      deallocate(X2)
      deallocate(X3)
      d_psi=matmul(A,X4)
      deallocate(X4)
    else
      d_psi=0
    endif
   else   ! ONLY 1 OBS
      oneobs=matmul(S,transpose(S))+R*dble(nrens)
      print *,'oneobs: ',oneobs(1,1)
      allocate (X3(nrobs,1))
      X3(1,1)=D(1,1)/oneobs(1,1)
     allocate (I_N(nrens,nrens))
      I_N=-1.0/dble(nrens)
      do i=1,nrens
         I_N(i,i)=I_N(i,i)+1.0
      enddo
      allocate(X4(nrens,1))
      X4=matmul(transpose(S),X3)
      deallocate(X3)
      allocate(X5(nrens,1))
      X5=matmul(I_N,X4)
      deallocate(X4)
      deallocate(I_N)
      d_psi=matmul(A,X5)
      deallocate(X5)
    end if

999 continue
    ! z0 uu vv [-2, 2]
    where (d_psi(1:2*nz+1,1) >  2.0) d_psi(1:2*nz+1,1) =  2.0
    where (d_psi(1:2*nz+1,1) < -2.0) d_psi(1:2*nz+1,1) = -2.0

    ! tt [-5, 5]
    where (d_psi(2*nz+2:3*nz+1,1) >  5.0) d_psi(2*nz+2:3*nz+1,1) =  5.0
    where (d_psi(2*nz+2:3*nz+1,1) < -5.0) d_psi(2*nz+2:3*nz+1,1) = -5.0

    ! ss [-6, 6]
    where (d_psi(3*nz+2:4*nz+1,1) >  6.0) d_psi(3*nz+2:4*nz+1,1) =  6.0
    where (d_psi(3*nz+2:4*nz+1,1) < -6.0) d_psi(3*nz+2:4*nz+1,1) = -6.0

!    if ((assim_mask == 0.0) .or. (T_surf < -1.5 .and. mlat > 40.) .or. &
!         (T_surf <  1.0 .and. mlat < -40.)) then
!      d_psi=0
!    end if

    psi=d_psi+psi

    psi4=REAL(psi(:,1))
    A4=REAL(A)
    R4=REAL(R)

end subroutine analysis2_letkf
