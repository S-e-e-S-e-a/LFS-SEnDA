module m_active_obs
contains
subroutine active_obs(i,j,obs,nrobs,lobs,wlobs,nobs,radius_sla,radius_sal,radius_tem,radius_sst,modlon,modlat)

! Calculate observations within range 'radius' of a model p-cell midpoint. 
! This is calculated for all points on grid.

  use mod_dimensions
  use mod_measurement
  use mod_angles
  use m_spherdist

  implicit none
  integer, intent(in)           :: i,j
  integer, intent(in)           :: nrobs
  real, intent(in)              :: radius_sla,radius_sal,radius_tem,radius_sst
  real, intent(in)              :: modlon(nx,ny),modlat(nx,ny)
  type(measurement), intent(in) :: obs(nrobs)
  logical, intent(out)          :: lobs(nrobs)
  real, intent(out)             :: wlobs(nrobs)
  integer, intent(out)          :: nobs

  real radius,z
  real lat0,lon0,pi_1,memdist
  integer m,ix,jy
  character(len=10) :: LocalFunc
  
  pi_1=4.*atan(1.)
  nobs=0
  do m=1,nrobs
     lon0 = ang180(obs(m)%lon+0.001) ! Add small number to avoid
     lat0 = obs(m)%lat+0.001         ! singularity in spherdist. Hrmph.
     memdist=spherdist(lon0,lat0,modlon(i,j),modlat(i,j))
      select case (trim(obs(m)%id))
             case ('TEM')
             radius=radius_tem
             case ('SAL')
             radius=radius_sal
             case ('SLA')
             radius=radius_sla
             case ('SST')
             radius=radius_sst
             case default
             print *,'No match in radius'
             print *,obs(m)%id
             stop
      end select 
    radius = radius*(cos(lat0*pi_1/180.)*0.6+0.4) !zyz change R
    LocalFunc='GC'
    if (LocalFunc == 'GC') then
    z = memdist / radius * 2.0
    
    if (z <= 1.) then
       nobs = nobs+1
       lobs(m)=.true.
       wlobs(m)=1. - 5./3.*z**2. + 5./8.*z**3. + 1./2.*z**4. - 1./4.*z**5. ! G-C Localization zyz
    elseif (z <= 2.) then
       nobs = nobs+1
       lobs(m)=.true.
       wlobs(m)=-2./3./z + 4. - 5.*z + 5./3.*z**2. + 5./8.*z**3. - 1./2.*z**4. + 1./12.*z**5.
    else
       lobs(m)=.false.
       wlobs(m)=0.
    endif
    elseif (LocalFunc == 'COS') then 
    z = memdist / radius
    if (z <= 1.) then
       nobs = nobs+1
       lobs(m)=.true.
       wlobs(m)=0.5*( 1.+cos(pi_1*z) )
    else
       lobs(m)=.false.
       wlobs(m)=0.
    endif
    else
      print *, 'LocalFunc Error: No this type!'
    endif
  enddo

end subroutine active_obs
end module m_active_obs
