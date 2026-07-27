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

  real radius
  real lat0,lon0,pi_1,memdist
  integer m,ix,jy

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

     if (spherdist(lon0,lat0,modlon(i,j),modlat(i,j)) <= radius) then
       nobs = nobs+1
       lobs(m)=.true.
!CY --- DEFINE WEIGHTING FUNCTION FOR ENSEMBLE OBS
!       wlobs(m)=1.         ! OLD choice : look out it is discontinuous 
!
!                          ! NEW choice : avoids discontinuity
!       wlobs(m)=0.5*( 1.+cos(pi_1*memdist/radius) )
!fww       wlobs(m)=(0.5*( 1.+cos(pi_1*memdist/radius) ))**2.
!lzt     
        wlobs(m)=(0.5*( 1.+cos(pi_1*memdist/radius) ))**0.5
!         wlobs(m)=1. 

!                          ! NEW choice : avoids discontinuity
!       if(memdist<=radius/2.) then
!       wlobs(m)=1.        
!       else
!       wlobs(m)=0.5*( 1.+cos(pi_1*(memdist-radius/2.)*2./radius) )
!       endif
!       obsindx(nobs) = m
     else
       lobs(m)=.false.
       wlobs(m)=0.
     endif
  enddo

end subroutine active_obs
end module m_active_obs
