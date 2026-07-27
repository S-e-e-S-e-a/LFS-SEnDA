module m_vert_loc
implicit none
contains
!===============================================================
! 
!   T_mod     : size=55
!   S_mod     : size=55
!   Lp        : scale fcator 0.5 kg/m3
!
!   Lz        : (ndim,nrobs)
!
!===============================================================

!---------------------------------------
! UNESCO 1983 seawater density
!---------------------------------------
real*8 function sw_rho(S, T)
    implicit none
    real*8, intent(in) :: S, T
    real*8 :: rho0, A, B

    rho0 = 999.842594d0 + 6.793952d-2*T - 9.095290d-3*T*T + 1.001685d-4*T**3  &
           - 1.120083d-6*T**4 + 6.536332d-9*T**5

    A = 8.24493d-1 - 4.0899d-3*T + 7.6438d-5*T*T - 8.2467d-7*T**3 + 5.3875d-9*T**4
    B = -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T*T

    sw_rho = rho0 + A*S + B*S**1.5d0 + 4.8314d-4*S*S
end function sw_rho

subroutine compute_Lz_density(nz, ndim, nrobs, T_mod, S_mod, dep_obs_ind, type_obs, Lz, Lc)
    implicit none

    integer, intent(in)    :: nz, nrobs, ndim
    real*8, intent(in)     :: T_mod(nz), S_mod(nz)
    integer, intent(in)    :: dep_obs_ind(nrobs)
    character(len=3), intent(in)    :: type_obs(nrobs)
    real*8, intent(out)    :: Lz(ndim, nrobs), Lc(ndim, nrobs)
    real*8  :: rho_mod(nz), Lz_single(nz,nrobs)
    real*8  :: drho, Lp
    integer :: i, j

    Lp=1.0
    Lc=1.0
    Lz=1.0
    Lz_single=1.0

 !   do i = 1, nz
 !     if ( S_mod(i) .gt. 9999. .or. T_mod(i) .gt. 9999. ) then
 !       rho_mod(i) = 9999.
 !     else
 !       rho_mod(i) = sw_rho(S_mod(i), T_mod(i))
 !     end if
 !   end do

 !   do j = 1, nrobs
 !       if ( trim(type_obs(j)) .ne. 'SLA') then
 !         do i = 1, nz
 !             drho    = abs(rho_mod(i) - rho_mod(dep_obs_ind(j)))
 !             Lz_single(i,j) = exp( - (drho / Lp)**2. )
 !         end do
 !       end if
 !   end do
 !   Lz(2:56,:)=Lz_single
 !   Lz(57:111,:)=Lz_single
 !   Lz(112:166,:)=Lz_single
 !   Lz(167:221,:)=Lz_single 
         
    do j = 1, nrobs
        if ( trim(type_obs(j)) == 'SLA' ) then
          do i=2,nz
          Lc(111+i,j) = .5
          Lc(166+i,j) = .5
          end do
       ! else if ( trim(type_obs(j)) == 'SST' .or. trim(type_obs(j)) == 'TEM' ) then   
       !   do i=1,nz
       !   Lc(166+i,j) = 0.2+0.8*(1-exp(-((rho_mod(i)-rho_mod(1))/Lp)**2))
       !   end do
       ! else if ( trim(type_obs(j)) == 'SAL') then
       !   do i=1,nz
       !   Lc(111+i,j) = 0.2+0.8*(1-exp(-((rho_mod(i)-rho_mod(1))/Lp)**2))
       !   end do
        end if
    end do

end subroutine compute_Lz_density

subroutine compute_Lz_depth(nz, ndim, nrobs, levels, dep_obs_ind, type_obs, Lz)
    implicit none

    integer, intent(in)    :: nz, nrobs, ndim
    real*4, intent(in)     :: levels(nz)
    integer, intent(in)    :: dep_obs_ind(nrobs)
    character(len=3), intent(in)    :: type_obs(nrobs)
    real*8, intent(out)    :: Lz(ndim,nrobs)
    
    real*8  :: Lz_single(nz, nrobs)
    real*8  :: dz, Lp, pi
    integer :: i, j

    pi = 4.*atan(1.)
    
    Lp=100.

    Lz=1.0
    Lz_single=1.0

    do j = 1, nrobs
       if ( trim(type_obs(j)) .ne. 'SLA') then
           do i = 1, nz
              dz    = abs(levels(i) - levels(dep_obs_ind(j)))
              Lz_single(i,j) = exp( - (dz / Lp)**2 )
           end do
        end if
    end do
    Lz(1,:)=Lz_single(1,:)
    Lz(2:nz+1,:)=Lz_single
    Lz(nz+2:2*nz+1,:)=Lz_single
    Lz(2*nz+2:3*nz+1,:)=Lz_single
    Lz(3*nz+2:4*nz+1,:)=Lz_single
 
end subroutine compute_Lz_depth

subroutine get_nz(type_obs, lon_obs, lat_obs, nobs, sqnz_obs)

implicit none
integer :: j, k, nobs
real*4  :: alpha
character(len=3), intent(in)  :: type_obs(nobs)
real*4, intent(in)            :: lon_obs(nobs), lat_obs(nobs)
real*4, intent(inout)         :: sqnz_obs(nobs)

sqnz_obs = 0.
alpha = 0.5
do j = 1, nobs
   if ( (trim(type_obs(j)) .ne. 'TEM') .and. (trim(type_obs(j)) .ne. 'SAL') ) then
     sqnz_obs(j) = 1.0
   else
     do k = 1, nobs
      if ( (trim(type_obs(j)) == trim(type_obs(k))) .and. (abs(lon_obs(j)-lon_obs(k)) < .1) .and. (abs(lat_obs(j)-lat_obs(k)) < .1) ) then
        sqnz_obs(j) = sqnz_obs(j) + 1.0
      endif
     enddo
   end if
   if ( sqnz_obs(j) > 1.001 ) then
     sqnz_obs(j) = (1+max(0.,sqnz_obs(j)-1.)*alpha)
   endif
enddo

end subroutine get_nz

end module m_vert_loc
